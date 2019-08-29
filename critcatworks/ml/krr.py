from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time, pathlib, sys
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import ase, ase.io
import logging, time
import numpy as np
from sklearn.kernel_ridge import KernelRidge
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, train_test_split

from critcatworks.database.extdb import update_machine_learning_collection
from critcatworks.database import read_descmatrix, write_descmatrix
import time

@explicit_serialize
class MLTask(FiretaskBase):
    """ 
    Machine Learning Task. It predicts the property of all uncomputed structures
    in the workflow. It is trained on all converged structures. Crossvalidation 
    is used to infer the optimal machine learning hyperparameters.
    Currently, only KRR (kernel ridge regression) is implemented. 

    A new document is added to the machine_learning collection. 
    
    Args:
        target_path (str) : absolute path to the target directory 
                            (needs to exist) on the computing resource.
    
    Returns:
        FWAction :  Firework action, updates fw_spec, possibly defuses
                    workflow upon failure.
    """
    _fw_name = 'MLTask'
    required_params = ['target_path']
    optional_params = []

    def run_task(self, fw_spec):
        target_path = self['target_path']

        METHOD = "krr"
        IS_PREDICT_FAILED = True
        N_CV = 5
        time_str = time.strftime("%Y-%m-%d-%H-%M")
        parent_folder_name = 'ml_' + METHOD + time_str
        parent_folder_path = target_path + "/" + parent_folder_name
        if not os.path.exists(parent_folder_path):
            os.makedirs(parent_folder_path)


        n_calcs_started = fw_spec["n_calcs_started"]
        calc_ids = fw_spec["temp"]["calc_ids"]
        is_converged_list = np.array(fw_spec["temp"]["is_converged_list"], dtype = 'bool')
        property_lst = fw_spec["temp"]["property"]
        descmatrix = read_descmatrix(fw_spec)
        workflow_id = fw_spec["workflow"]["_id"]
        workflow_parameters = fw_spec["workflow"]["parameters"]

        enumerated_ids = np.arange(len(calc_ids))

        ### PREPARE, CHECK ###
        is_converged_ids = np.array(calc_ids)[is_converged_list]

        print("calc_ids", calc_ids)
        print("is_converged_list", is_converged_list)
        print("is_converged_ids", is_converged_ids)

        finished_calc_ids = np.array(calc_ids)[:n_calcs_started]
        simulation_ids_training,  training_ids, _ = np.intersect1d(finished_calc_ids, is_converged_ids, return_indices = True)

        print("training_ids", training_ids)

        # magic number due to n-fold CV
        if training_ids.shape[0] < N_CV:
            logging.warning('Exiting workflow Problem detected: ' + 
                'Too few datapoints to learn from! Something might be wrong with your DFT calculations!')
            return FWAction(defuse_children=True, update_spec = fw_spec)

        if IS_PREDICT_FAILED == True:
            simulation_ids_predict = np.array(calc_ids)[is_converged_list == 0]
            to_predict_ids = np.array(enumerated_ids)[is_converged_list == 0]
        else:
            to_predict_ids = np.array(calc_ids)[n_calcs_started:]

        
        features = descmatrix[training_ids]
        to_predict_features = descmatrix[to_predict_ids]
        labels = np.array(property_lst)[training_ids]

        print("features", features.shape)
        print("features", labels, labels.shape)
        print("training_ids", training_ids)
        print("to_predict_ids", to_predict_ids)
        ### RUN ###

        if METHOD == "krr":
            ml_results = ml_krr(features, labels, training_ids, 
                to_predict_features, to_predict_ids, 
                is_scaled = False, n_cv = N_CV, path = parent_folder_path)
        else:
            logging.warning('Exiting workflow Problem detected: ' + 
                'machine learning method  ' + METHOD + '  not implemented')
            return FWAction(defuse_children=True)

        ### DATABASE ###
        # update machine learning data
        # translate to simulation ids
        
        ml_results["ids_train"] = simulation_ids_training[ml_results["ids_train"]]
        ml_results["ids_test"] = simulation_ids_training[ml_results["ids_test"]]
        ml_results["ids_predicted"] = simulation_ids_predict

        # update external database
        dct = update_machine_learning_collection(METHOD, extdb_connect = fw_spec["extdb_connect"], workflow_id = workflow_id, 
            method_params = ml_results["method_params"], 
            descriptor = workflow_parameters["descriptor"],
            descriptor_params = workflow_parameters["descriptor_params"],
            training_set = ml_results["ids_train"].tolist(), 
            validation_set = [], 
            test_set = ml_results["ids_test"].tolist(), 
            prediction_set = simulation_ids_predict.tolist(),
            metrics_training = ml_results["metrics_training"], 
            metrics_validation = ml_results["metrics_validation"],
            metrics_test = ml_results["metrics_test"],
            output = ml_results["output"]
            )
        
        # update internal workflow data
        update_spec = fw_spec
        machine_learning_id = dct["_id"]
        update_spec["machine_learning"][str(machine_learning_id)] = dct

        # update temp workflow data
        update_spec["temp"]["last_machine_learning_id"] = machine_learning_id
        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)


def get_mae(target_path):
    """ 
    Creates Firework from MLTask. It predicts the property of all uncomputed structures
    in the workflow. It is trained on all converged structures. Crossvalidation 
    is used to infer the optimal machine learning hyperparameters.
    Currently, only KRR (kernel ridge regression) is implemented. 

    A new document is added to the machine_learning collection. 
    
    Args:
        target_path (str) : absolute path to the target directory 
                            (needs to exist) on the computing resource.
    
    Returns:
        FWAction :  Firework action, updates fw_spec, possibly defuses
                    workflow upon failure.
    """
    firetask1  = MLTask(target_path=target_path)
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'MLTask'},
             name = 'MLWork')
    return fw


def ml_krr(features, labels, train_test_ids, to_predict_features, to_predict_ids, 
        alpha_list= np.logspace(-1, -9, 9),
        gamma_list = np.logspace(-1, -9, 9), 
        kernel_list = ['rbf'], 
        sample_size=0.8,
        is_scaled = False,
        n_cv = 5,
        path = "."):
    """
    Helper function to estimate the generalization error (MAE, MSE). The hyperparameters alpha and gamma are
    by default scanned on a logarithmic scale. The data set is split randomly into training and test set.
    The ratio of the split is defined by sample_size.
    The training set is used for cross validation.

    Args:
        features (2D ndarray) : descriptor input for the machine learning algorithm for training/testing
        labels (1D ndarray) :   property labels for the machine learning algorithm for training/testing
        train_test_ids (1D ndarray) :   pythonic ids (of features and labels) for training and
                                        testing. 
        to_predict_features (1D ndarray) :  descriptor input for the machine learning algorithm 
                                            for prediction
        to_predict_ids (1D ndarray) :   pythonic ids (of features and labels) ommited from training and
                                        testing. 
        alpha_list (lsit) :     Regularization parameter. Defaults to np.logspace(-1, -9, 9)
        gamma_list (list) :     Kernel function scaling parameter. Defaults to np.logspace(-1, -9, 9)
        kernel_list (list) :    List of kernel functions (see sklearn documentation for options).
                                Defaults to ['rbf']
        sample_size (float) : The ratio of the training-test split is defined by this. Defaults to 0.8
        is_scaled (bool) : If set to True, the features are scaled. Defaults to False
        n_cv (int) :    Number of cross-validation splits. Defaults to 5
        path (str) :    path whereto to write the machine learning output. Defaults to the
                        current working directory

    Returns:
        dict :  machine learning results with the following keys:
                ids_train, ids_test, ids_predicted, method_params, 
                output (.label_predicted, .label_train, .label_test),
                metrics_test, metrics_validation, metrics_training
    """
    # load, split and scale data
    x_train, x_test, y_train, y_test, ids_train, ids_test = split_scale_data(features, labels, train_test_ids, sample_size, is_scaled)

    # Create kernel linear ridge regression object
    learner = GridSearchCV(KernelRidge(kernel='rbf'), n_jobs = 8, cv=n_cv,
                param_grid={"alpha": alpha_list, "gamma": gamma_list, 
                "kernel": kernel_list}, scoring = 'neg_mean_absolute_error', return_train_score=True)

    t_ml0 = time.time()
    learner.fit(x_train, y_train)
    t_ml1 = time.time()
    print("ml time", str(t_ml1 - t_ml0))

    # getting best parameters
    learner_best = learner.best_estimator_

    mae, mse, y_pred, train_y_pred, learner_best = predict_and_error(learner_best, x_test, x_train, y_test)

    # predict remaining datapoints
    y_to_predict = learner_best.predict(to_predict_features)

    ### OUTPUT ###
    write_output(learner, sample_size, "krr", mae, mse, "param", 
        ids_test, y_test, y_pred, 
        ids_train, y_train, train_y_pred,
        to_predict_ids, y_to_predict,
        path,
        )

    ml_results = {
        "ids_train" : ids_train,
        "ids_test" : ids_test,
        "ids_predicted" : to_predict_ids,
        "method_params" : learner.best_params_,
        "output" : {"label_predicted" : y_to_predict.tolist(), "label_train" : train_y_pred.tolist(), "label_test" : y_pred.tolist()},
        "metrics_test" : {"mae" : mae, "mse" : mse},
        "metrics_validation" : {"mae" : -1 * learner.cv_results_['mean_test_score'][learner.best_index_],
            "std" : learner.cv_results_['std_test_score'][learner.best_index_]},
        "metrics_training" : {"mae" : -1 * learner.cv_results_['mean_train_score'][learner.best_index_],
            "std" : learner.cv_results_['std_train_score'][learner.best_index_]},
        }


    return ml_results

def scale_data(x_train, x_test, is_mean=True):
    """
    Helper function to the scale the data with respect
    to the mean.

    Args:
        x_train (2D ndarray) : training data which is scaled
        x_test (2D ndarray) : test data scaled accordingly
        is_mean (bool) :    if set to False, scaled between 0 and 1.
                            Otherwise, scaled centered around the mean of x_train.

    Returns:
        tuple : the scaled arrays x_train, x_test
    """
    # Scale
    scaler = StandardScaler(with_mean=is_mean)  
        # fit only on training data
    scaler.fit(x_train)  
    x_train = scaler.transform(x_train)  
        # apply same transformation to test data
    x_test = scaler.transform(x_test)  
    return x_train, x_test

def split_scale_data(x_data, y_data, ids_data, sample_size, is_scaled):
    """
    Helper function to split and scale the data.

    Args:
        x_data (2D ndarray) : features of the training and test data
        y_data (1D ndarray) : labels of the training and test data
        ids_data (1D ndarray) : complete list of ids of datapoints used for training and testing
        sample_size (float) : The ratio of the training-test split is defined by this
        is_scaled (bool) : True scales the features centered around the mean

    Returns:
        tuple : ndarrays    x_train, x_test (split from x_data) 
                            y_train, y_test (split from y_data)
                            ids_train, ids_test (split from ids_data)
    """
    x_train, x_test, y_train, y_test, ids_train, ids_test = train_test_split(x_data, y_data, ids_data, test_size = 1 - sample_size)

        # scale
    if is_scaled:
        x_train, x_test = scale_data(x_train, x_test, is_mean=True)
    return x_train, x_test, y_train, y_test, ids_train, ids_test

def predict_and_error(learner, x_test, x_train, y_test):
    """
    Helper function to predict the property on a training and a test set.


    Args:
        learner (sklearn.learner) : learner object with which the training set was fitted
        x_test (2D ndarray) :   test data scaled accordingly
        x_train (2D ndarray) :  training data which is scaled
        y_test (1D ndarray) :   labels of the test data

    Returns:
        tuple : mae, mse, y_pred, train_y_pred, learner
    """
    y_pred = learner.predict(x_test)

    # run also on training set
    train_y_pred = learner.predict(x_train)

    # errors
    mae = np.absolute(y_pred - y_test)
    mse = mae ** 2

    mae = np.mean(mae)
    mse = np.mean(mse)
    return mae, mse, y_pred, train_y_pred, learner


def write_output(learner, sample_size, ml_method, mae, mse, runtype, 
    ids_test, y_test, y_pred, ids_train, y_train, train_y_pred,
    to_predict_ids, y_to_predict, path):
    """
    Helper function to write the machine learning output.
    """
    ### OUTPUT ###
    # y_test vs y_predict
    y_tmp = np.array([ids_test, y_test, y_pred])
    y_compare = np.transpose(y_tmp)

    np.savetxt(path + "/" + ml_method + str("_") + runtype + "_size" + str(sample_size) + ".predictions", y_compare, 
        header = "###ids_test   y_test    y_pred")

    # also y_train vs. y_pred_train
    y_tmp = np.array([ids_train, y_train, train_y_pred])
    y_compare = np.transpose(y_tmp)

    np.savetxt(path + "/" + ml_method + str("_") + runtype + "_size" + str(sample_size) + ".trainset_predictions", y_compare, 
        header = "###ids_train   y_train    train_y_pred")

    # also to_predict_ids and y_to_predict
    y_tmp = np.array([to_predict_ids, y_to_predict])
    y_compare = np.transpose(y_tmp)

    np.savetxt(path + "/" + ml_method + str("_") + runtype + "_size" + str(sample_size) + ".remaining_predictions", y_compare, 
        header = "###ids_remaining   y_pred_remaining")

    with open(path + "/" + ml_method + str("_") + runtype + "_size" + str(sample_size) + ".out", "w") as f:
        f.write("MAE " + str(mae) + "\n")
        f.write("MSE " + str(mse) + "\n")
        f.write("\n")
        if runtype == "param":
            f.write("Best parameters of " + ml_method + ": \n")
            f.write(str(learner.best_params_) + "\n")
            f.write("Errors of best parameters: \n")



            f.write("Grid scores on validation set:" + "\n")
            f.write("\n")
            means = learner.cv_results_['mean_test_score']
            stds = learner.cv_results_['std_test_score']
            for mean, std, params in zip(means, stds, learner.cv_results_['params']):
                f.write("%0.4f (+/-%0.04f) for %r"
                      % (mean, std * 2, params) + "\n")
            f.write("\n")
            f.write("Grid scores on train set:" + "\n")
            f.write("\n")        
            means = learner.cv_results_['mean_train_score']
            stds = learner.cv_results_['std_train_score']        
            for mean, std, params in zip(means, stds, learner.cv_results_['params']):
                f.write("%0.4f (+/-%0.04f) for %r"
                      % (mean, std * 2, params) + "\n")
    return None
