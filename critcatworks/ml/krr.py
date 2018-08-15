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
import sklearn


@explicit_serialize
class MLTask(FiretaskBase):
    """ 
    Task to update database from converged chunk
    of calculations.
    """

    _fw_name = 'MLTask'
    required_params = []
    optional_params = []

    def run_task(self, fw_spec):

        logging.info("ML not implemented yet")

        n_calcs_started = fw_spec["n_calcs_started"]
        ranked_ids = fw_spec["fps_ranking"][:n_calcs_started]
        to_predict_ids = fw_spec["fps_ranking"][n_calcs_started:]


        descmatrix = np.array(fw_spec["descmatrix"])
        toten = np.array(fw_spec["adsorbate_energies_list"])
        
        features = descmatrix[ranked_ids]
        to_predict_features = descmatrix[ranked_ids]
        labels = toten[ranked_ids]

        mae, y_to_predict, krr_parameters = ml_krr(features, labels, ranked_ids, to_predict_features, to_predict_ids, is_scaled = False)

        update_spec = fw_spec
        update_spec["mae"] = mae
        update_spec["best_krr_parameters"] = krr_parameters
        update_spec["ids_predicted"] = to_predict_ids
        update_spec["predicted_energies"] = y_to_predict

        return FWAction(update_spec=update_spec)


def get_mae():
    firetask1  = MLTask()
    fw = Firework([firetask1])
    return fw


def ml_krr(features, labels, train_test_ids, to_predict_features, to_predict_ids, 
        alpha_list= np.logspace(-1, -9, 9),
        gamma_list = np.logspace(-1, -9, 9), 
        kernel_list = ['rbf'], 
        sample_size=0.8,
        is_scaled = False):
    
    # load, split and scale data
    x_train, x_test, y_train, y_test, ids_train, ids_test = split_scale_data(x_data, y_data, ids, sample_size, is_scaled)

    # Create kernel linear ridge regression object
    learner = GridSearchCV(KernelRidge(kernel='rbf'), n_jobs = 8, cv=5,
                param_grid={"alpha": alpha_list, "gamma": gamma_list, 
                "kernel": kernel_list}, scoring = 'neg_mean_absolute_error')

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
        to_predict_ids, y_to_predict
        )

    return mae, y_to_predict, learner.best_params_


def scale_data(x_train, x_test, is_mean=True):
    # Scale
    scaler = StandardScaler(with_mean=is_mean)  
        # fit only on training data
    scaler.fit(x_train)  
    x_train = scaler.transform(x_train)  
        # apply same transformation to test data
    x_test = scaler.transform(x_test)  
    return x_train, x_test

def split_scale_data(x_datafile, y_datafile, sample_size, is_scaled):
    
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
    to_predict_ids, y_to_predict):
    ### OUTPUT ###
    # y_test vs y_predict
    y_tmp = np.array([ids_test, y_test, y_pred])
    y_compare = np.transpose(y_tmp)

    np.savetxt(ml_method + str("_") + runtype + "_size" + str(sample_size) + ".predictions", y_compare, 
        header = "###ids_test   y_test    y_pred")

    # also y_train vs. y_pred_train
    y_tmp = np.array([ids_train, y_train, train_y_pred])
    y_compare = np.transpose(y_tmp)

    np.savetxt(ml_method + str("_") + runtype + "_size" + str(sample_size) + ".trainset_predictions", y_compare, 
        header = "###ids_train   y_train    train_y_pred")

    # also to_predict_ids and y_to_predict
    y_tmp = np.array([to_predict_ids, y_to_predict])
    y_compare = np.transpose(y_tmp)

    np.savetxt(ml_method + str("_") + runtype + "_size" + str(sample_size) + ".remaining_predictions", y_compare, 
        header = "###ids_remaining   y_pred_remaining")



    with open(ml_method + str("_") + runtype + "_size" + str(sample_size) + ".out", "w") as f:
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