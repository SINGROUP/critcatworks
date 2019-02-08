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

@explicit_serialize
class MLTask(FiretaskBase):
    """ 
    Task to update database from converged chunk
    of calculations.
    """

    _fw_name = 'MLTask'
    required_params = ['target_path']
    optional_params = []

    def run_task(self, fw_spec):
        target_path = self['target_path']

        METHOD = "krr"
        IS_PREDICT_FAILED = True
        N_CV = 5
        parent_folder_name = 'ml_' + METHOD
        parent_folder_path = target_path + "/" + parent_folder_name
        if not os.path.exists(parent_folder_path):
            os.makedirs(parent_folder_path)


        n_calcs_started = fw_spec["n_calcs_started"]
        calc_ids = fw_spec["temp"]["calc_ids"]
        is_converged_list = np.array(fw_spec["temp"]["is_converged_list"], dtype = 'bool')
        property_lst = fw_spec["temp"]["property"]
        descmatrix = np.array(fw_spec["temp"]["descmatrix"])
        workflow_id = fw_spec["workflow"]["_id"]
        workflow_parameters = fw_spec["workflow"]["parameters"]

        enumerated_ids = np.arange(len(calc_ids))

        ### PREPARE, CHECK ###
        is_converged_ids = np.array(calc_ids)[is_converged_list]

        print("calc_ids", calc_ids)
        print("is_converged_list", is_converged_list)
        print("is_converged_ids", is_converged_ids)

        finished_calc_ids = np.array(calc_ids)[:n_calcs_started]
        simulation_ids_training, _ , training_ids = np.intersect1d(finished_calc_ids, is_converged_ids, return_indices = True)
        #training_ids = np.array(enumerated_ids)[np.intersect1d(finished_calc_ids, is_converged_ids)]

        print("training_ids", training_ids)

        # magic number due to n-fold CV
        if training_ids.shape[0] < N_CV:
            logging.warning('Exiting workflow Problem detected: ' + 
                'Too few datapoints to learn from! Something might be wrong with your DFT calculations!')
            return FWAction(defuse_workflow=True)

        if IS_PREDICT_FAILED == True:
            simulation_ids_predict = np.array(calc_ids)[1 - is_converged_list]
            to_predict_ids = np.array(enumerated_ids)[1 - is_converged_list]
        else:
            to_predict_ids = np.array(calc_ids)[n_calcs_started:]

        
        features = descmatrix[training_ids]
        to_predict_features = descmatrix[to_predict_ids]
        labels = np.array(property_lst)[training_ids]

        ### RUN ###

        if METHOD == "krr":
            #mae, mse, y_to_predict, method_params = ml_krr(features, labels, training_ids, 
            ml_results = ml_krr(features, labels, training_ids, 
                to_predict_features, to_predict_ids, 
                is_scaled = False, n_cv = N_CV, path = parent_folder_path)
        else:
            logging.warning('Exiting workflow Problem detected: ' + 
                'machine learning method  ' + METHOD + '  not implemented')
            return FWAction(defuse_workflow=True)

        ### DATABASE ###
        # update machine learning data
        # translate to simulation ids
        
        ml_results["ids_train"] = simulation_ids_training[ml_results["ids_train"]]
        ml_results["ids_test"] = simulation_ids_training[ml_results["ids_test"]]
        ml_results["ids_predicted"] = simulation_ids_predict

        # update external database
        dct = update_machine_learning_collection(METHOD, workflow_id = workflow_id, 
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

        #update_spec["mae"] = mae
        #update_spec["best_krr_parameters"] = krr_parameters
        #update_spec["ids_predicted"] = to_predict_ids
        #update_spec["predicted_energies"] = y_to_predict

        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)


def get_mae(target_path):
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

    #return mae, mse, y_to_predict, learner.best_params_
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
    # Scale
    scaler = StandardScaler(with_mean=is_mean)  
        # fit only on training data
    scaler.fit(x_train)  
    x_train = scaler.transform(x_train)  
        # apply same transformation to test data
    x_test = scaler.transform(x_test)  
    return x_train, x_test

def split_scale_data(x_data, y_data, ids_data, sample_size, is_scaled):
    
    x_train, x_test, y_train, y_test, ids_train, ids_test = train_test_split(x_data, y_data, ids_data, test_size = 1 - sample_size)

        # scale
    if is_scaled:
        x_train, x_test = scale_data(x_train, x_test, is_mean=True)
    return x_train, x_test, y_train, y_test, ids_train, ids_test

def predict_and_error(learner, x_test, x_train, y_test):

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
