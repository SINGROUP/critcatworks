# utility script to do post-run  machine-learning analysis
# runs ml multiple times in order to get an estimate of
# the variance
# requires workflow_id of the mongodb database
# default is the critcat database, the flag --is-test redirects to the test database

from critcatworks.database.format import atoms_dict_to_ase, ase_to_atoms_dict
from critcatworks.database.extdb import fetch_simulations
import ase, ase.io
from ase.visualize import view
import pymongo
from pprint import pprint as pp
import argparse, json, time
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
from get_ml_from_workflow import get_external_database
from get_ml_from_workflow import get_ml
from get_ml_from_workflow import get_property
from critcatworks.ml.krr import split_scale_data, predict_and_error, write_output, ml_krr
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.preprocessing import StandardScaler

import os

def krr_multirun_task(simulation_ids_training, descmatrix, property_lst):
    METHOD = "krr"
    IS_PREDICT_FAILED = True
    N_CV = 5
    time_str = time.strftime("%Y-%m-%d-%H-%M")
    parent_folder_path = "krr_multirun"
    if not os.path.exists(parent_folder_path):
        os.makedirs(parent_folder_path)

    calc_ids = np.array(CALC_IDS)
    print("calc ids" , calc_ids.shape, calc_ids[0])
    simulation_ids_training = np.array(simulation_ids_training)
    simulation_ids_predict = np.setdiff1d(calc_ids, simulation_ids_training)

    enumerated_ids = np.arange(len(calc_ids))
    print("enum ids", enumerated_ids.shape)

    ### PREPARE, CHECK ###
    print("simulation_ids_training", simulation_ids_training.shape, simulation_ids_training[0])
    simulation_ids_training,  training_ids, _ = np.intersect1d(calc_ids, simulation_ids_training, return_indices = True)

    to_predict_ids = np.setdiff1d(enumerated_ids, training_ids)

    features = descmatrix[training_ids]
    to_predict_features = descmatrix[to_predict_ids]
    labels = np.array(property_lst)[training_ids]

    print("features", features.shape)
    print("labels", labels.shape)
    print("training_ids", training_ids.shape)
    print("to_predict_ids", to_predict_ids.shape)
    ### RUN ###

    RUNS = 10
    N_CV = 5
    mae_test_lst = []
    mae_validation_lst = []
    mae_training_lst = []
    for run in range(RUNS):
        print("ML run number", run)

        ml_results = ml_krr(features, labels, training_ids, 
            to_predict_features, to_predict_ids, 
            is_scaled = False, n_cv = N_CV, path = parent_folder_path)
        # store mae test 
        mae_test_lst.append(ml_results["metrics_test"]["mae"])
        mae_validation_lst.append(ml_results["metrics_validation"]["mae"])
        mae_training_lst.append(ml_results["metrics_training"]["mae"])
        method_params = ml_results["method_params"]

    mae_test_arr = np.array(mae_test_lst)
    mae_validation_arr = np.array(mae_validation_lst)
    mae_training_arr = np.array(mae_training_lst)
    print("test multiple MAE", mae_test_arr)
    averaged_test_mae = np.mean(mae_test_arr)
    std_test_mae = np.std(mae_test_arr)
    averaged_validation_mae = np.mean(mae_validation_arr)
    std_validation_mae = np.std(mae_validation_arr)
    averaged_training_mae = np.mean(mae_training_arr)
    std_training_mae = np.std(mae_training_arr)
    
    # once run full set with optimized parameters
    method_params
    y_to_predict, train_y_predict = krr_optimized_hyperparameters(features, labels, training_ids,
        to_predict_features, to_predict_ids,
        method_params, is_scaled = False, path = parent_folder_path)

    ### DATABASE ###
    # update machine learning data
    # translate to simulation ids
    
    ml_results["ids_train"] = np.array(calc_ids)[ml_results["ids_train"]]
    ml_results["ids_test"] = np.array(calc_ids)[ml_results["ids_test"]]
    ml_results["ids_predicted"] = simulation_ids_predict

    dct = dict(
        method = METHOD, 
        #workflow_id = workflow_id, 
        method_params = ml_results["method_params"], 
        #descriptor = workflow_parameters["descriptor"],
        #descriptor_params = workflow_parameters["descriptor_params"],
        training_set = simulation_ids_training.tolist(), 
        validation_set = [],
        test_set = [],
        prediction_set = simulation_ids_predict.tolist(),
        metrics_training = {"mae" : averaged_training_mae, "std" : std_training_mae},
        metrics_validation ={"mae" : averaged_validation_mae, "std" : std_validation_mae},
        metrics_test = {"mae" : averaged_test_mae, "std" : std_test_mae},
        output = {"labels_predict" : y_to_predict.tolist(), "labels_train" : train_y_predict},
        )
    return dct  # machine learning dict, overwrite info later 


def krr_optimized_hyperparameters(features, labels, training_ids, to_predict_features, to_predict_ids, method_params, is_scaled = False, path = "."):
    sample_size = 1.0
    # load, split and scale data
    #x_train, x_test, y_train, y_test, ids_train, ids_test = split_scale_data(features, labels, train_test_ids, sample_size, is_scaled)

    # Create kernel linear ridge regression object
    #learner = GridSearchCV(KernelRidge(kernel='rbf'), n_jobs = 8, cv=n_cv,
    #            param_grid={"alpha": alpha_list, "gamma": gamma_list, 
    #            "kernel": kernel_list}, scoring = 'neg_mean_absolute_error', return_train_score=True)
    learner = KernelRidge(kernel = method_params["kernel"], alpha = method_params["alpha"], gamma = method_params["gamma"])

    t_ml0 = time.time()
    learner.fit(features, labels)
    t_ml1 = time.time()
    print("ml time", str(t_ml1 - t_ml0))

    # getting best parameters
    #learner_best = learner.best_estimator_

    #mae, mse, y_pred, train_y_pred, learner_best = predict_and_error(learner_best, x_test, x_train, y_test)
    train_y_predict = learner.predict(features)

    # predict remaining datapoints
    y_to_predict = learner.predict(to_predict_features)

    ### OUTPUT ###
    # modify writing output
    ids_test = []
    y_test = []
    y_pred = []

    write_output(learner, sample_size, "krr", None, None, "run",
        ids_test, y_test, y_pred,
        training_ids, labels, train_y_predict,
        to_predict_ids, y_to_predict,
        path,
        )
    return y_to_predict, train_y_predict


def krr_fps_ordered_run(simulation_ids_training, simulation_ids_test, descmatrix, property_lst):
    METHOD = "krr"
    IS_PREDICT_FAILED = True
    N_CV = 5
    time_str = time.strftime("%Y-%m-%d-%H-%M")
    parent_folder_path = "krr_multirun"
    if not os.path.exists(parent_folder_path):
        os.makedirs(parent_folder_path)

    calc_ids = np.array(CALC_IDS)
    print("calc ids" , calc_ids.shape, calc_ids[0])
    simulation_ids_training = np.array(simulation_ids_training)
    simulation_ids_training = np.array(simulation_ids_test)
    simulation_ids_predict = np.setdiff1d(calc_ids, simulation_ids_training)

    enumerated_ids = np.arange(len(calc_ids))
    print("enum ids", enumerated_ids.shape)

    ### PREPARE, CHECK ###
    print("simulation_ids_training", simulation_ids_training.shape, simulation_ids_training[0])
    simulation_ids_training,  training_ids, _ = np.intersect1d(calc_ids, simulation_ids_training, return_indices = True)
    simulation_ids_test,  test_ids, _ = np.intersect1d(calc_ids, simulation_ids_test, return_indices = True)

    to_predict_ids = np.setdiff1d(enumerated_ids, training_ids)

    features = descmatrix[training_ids]
    to_predict_features = descmatrix[to_predict_ids]
    labels = np.array(property_lst)[training_ids]
    features_test = descmatrix[test_ids]
    labels_test = np.array(property_lst)[test_ids]

    print("features", features.shape)
    print("labels", labels.shape)
    print("training_ids", training_ids.shape)
    print("to_predict_ids", to_predict_ids.shape)
    ### RUN ###

    N_CV = 5

    #ml_results = ml_krr(features, labels, training_ids, 
    #    to_predict_features, to_predict_ids, 
    #    is_scaled = False, n_cv = N_CV, path = parent_folder_path)

    # Create kernel linear ridge regression object
    alpha_list= np.logspace(-1, -9, 9) 
    gamma_list = np.logspace(-1, -9, 9) 

    learner = GridSearchCV(KernelRidge(kernel='rbf'), n_jobs = 8, cv=N_CV,
                param_grid={"alpha": alpha_list, "gamma": gamma_list}, 
                scoring = 'neg_mean_absolute_error', return_train_score=True)

    t_ml0 = time.time()
    learner.fit(features, labels)
    t_ml1 = time.time()
    print("ml time", str(t_ml1 - t_ml0))

    # getting best parameters
    learner_best = learner.best_estimator_

    mae, mse, y_pred, train_y_pred, learner_best = predict_and_error(learner_best, features_test, features, labels_test)
    train_y_predict = learner_best.predict(features)

    # predict remaining datapoints
    y_to_predict = learner_best.predict(to_predict_features)

    ### OUTPUT ###
    # modify writing output

    write_output(learner, 1.0, "krr_fps_ordered", None, None, "param",
        test_ids, labels_test, y_pred,
        training_ids, labels, train_y_predict,
        to_predict_ids, y_to_predict,
        parent_folder_path,
        )
    return mae


def get_learning_curve(workflow_id, is_test, descmatrix, property_lst):
    # train, validation and test MAE
    training_sizes = []
    mae_train = []
    mae_validation = []
    mae_test = []
    std_train = []
    std_validation = []
    std_test = []

    mae_fps_ordered = []

    if True:
        ml = get_ml(workflow_id, is_test)
        count = 0
        for entry in ml:

            # store previous set
            full_set = []
            full_set.extend(entry["training_set"])
            full_set.extend(entry["test_set"])


            ml_results = krr_multirun_task(np.array(full_set), descmatrix, property_lst)

            train_size = len(ml_results["training_set"])
            training_sizes.append(train_size)
            test_size = len(ml_results["test_set"])
            print("train size", train_size, "test_size", test_size)
            mae_train.append(ml_results["metrics_training"]["mae"] * 27.2114)
            mae_validation.append(ml_results["metrics_validation"]["mae"] * 27.2114)
            mae_test.append(ml_results["metrics_test"]["mae"] * 27.2114)

            std_train.append(ml_results["metrics_training"]["std"] * 27.2114)
            std_validation.append(ml_results["metrics_validation"]["std"] * 27.2114)
            std_test.append(ml_results["metrics_test"]["std"] * 27.2114)


            if count == 0:
                count += 1
                previous_set = full_set.copy()
                continue
            else:
                # mae test future
                # fps orderered training and test set
                simulation_ids_training, _ , _ = np.intersect1d(full_set, previous_set, return_indices = True)
                simulation_ids_test = np.setdiff1d(full_set, previous_set)

                previous_set = full_set.copy()
                mae = krr_fps_ordered_run(simulation_ids_training, simulation_ids_test, descmatrix, property_lst)
                mae_fps_ordered.append(mae * 27.2114)
                print("#########################")
                print("MAE fps ordered", mae * 27.2114)
                print("#########################")


    def draw_func():
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        #errorbar(x, y, e, linestyle='None', marker='^')
        #ax.plot(training_sizes, mae_train, label = 'training set')
        #ax.plot(training_sizes, mae_validation, label = 'validation set')
        #ax.plot(training_sizes, mae_test, label = 'test set')
        ax.plot(training_sizes[:-1], mae_fps_ordered, label = 'test set fps ordered')
        ax.errorbar(training_sizes, mae_train, std_train, label = 'training set')
        ax.errorbar(training_sizes, mae_validation, std_validation, label = 'validation set')
        ax.errorbar(training_sizes, mae_test, std_test, label = 'test set')
        ax.xaxis.label.set_fontsize(20)
        ax.set_xlabel('training set size', size=20)
        #ax.set_ylabel(r'$\alpha^2$ F [Ha]', size=20)
        ax.set_ylabel('MAE [eV]', size=20)
        ax.ticklabel_format(style='sci', axis='y', size =20, scilimits=(0,0))
        ax.legend(loc = 'upper right')
        plt.savefig('learning_curve_multirun.png')
        plt.show()
        return None

    y_tmp = np.array([training_sizes, mae_train, std_train, mae_validation, std_validation, mae_test, std_test])
    y_compare = np.transpose(y_tmp)
    np.savetxt("summary_multirun_krr.dat", y_compare,
        header = "###training_sizes, mae_train, std_train, mae_validation, std_validation, mae_test, std_test")
    y_tmp = np.array([training_sizes[:-1], mae_fps_ordered])
    y_compare = np.transpose(y_tmp)
    np.savetxt("summary_multirun_krr.fps_ordered_test.dat", y_compare,
        header = "###training_sizes[:-1], mae_fps_ordered")

    draw_func()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('ids', metavar='N', type=int, nargs='+',
                    help='workflow id')
    parser.add_argument('--test', action='store_true')


    # read FW.json
    with open ("FW.json", 'r') as f:
       js = json.load(f)

    fw_spec = js["spec"]
    CALC_IDS = fw_spec["temp"]["calc_ids"]
    PROPERTY_LST = fw_spec["temp"]["property"]


    # descmatrix needs to be loaded beforehand
    DESCMATRIX = np.load("descmatrix_2019-07-27-14-52.npy")

    args = parser.parse_args()
    print(args.ids)
    is_test = args.test
    for idx in args.ids:
        ml = get_ml(idx, is_test)
        for entry in ml:
            train_size = len(entry["training_set"])
            test_size = len(entry["test_set"])
            print("id", entry["_id"])
            print("train size", train_size, "test_size", test_size)
            print("method parameters", entry["method_params"])
            #entry["training_set"] = 0
            #entry["test_set"] = 0
            #entry["output"] = 0
            #pp(entry)


    # parity plot of last entry
    #get_parity_plot(entry)


    get_learning_curve(args.ids[0], is_test, DESCMATRIX, PROPERTY_LST)

