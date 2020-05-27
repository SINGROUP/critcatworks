# get machine learning information from a
# given workflow_id
# parity plot of the last machine learning entry is plotted
# finally, a learning curve is plotted

from critcatworks.database.format import atoms_dict_to_ase, ase_to_atoms_dict
from critcatworks.database.extdb import fetch_simulations
import ase, ase.io
from ase.visualize import view
import pymongo
from pprint import pprint as pp
import argparse
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib

def get_external_database(**extdb_connect):
    #extdb_connect["username"] = username
    #extdb_connect["password"] = password
    extdb_connect["host"] = extdb_connect.get("host",
        "nanolayers.dyndns.org:27017")
    extdb_connect["db_name"] = extdb_connect.get("db_name",
        "testdb")        
    extdb_connect["authsource"] = extdb_connect.get("authsource",
        extdb_connect["db_name"])
    CLIENT = pymongo.MongoClient(extdb_connect["host"], 
            username = extdb_connect["username"],
            password = extdb_connect["password"],
            authSource = extdb_connect["authsource"])
    db = CLIENT[extdb_connect["db_name"]]
    return db, extdb_connect

def get_simulations_from_elements(el1, el2, workflow_id = 45, is_test = False):
    if is_test:
        db_name = "testdb"
    else:
        db_name = "ncdb"

    ext_db = get_external_database(
        username = "myusername",
        password = "mypassword",
        db_name = db_name,
        )
    # subset tests for metal elements plus common adsorbate light elements
    cursor = ext_db.simulations.find(
        {'workflow_id' : workflow_id, 
        "$expr":{"$setIsSubset":["$atoms.numbers",[el1,el2, 0,1,6,7,8,9]]}
        }, 
        {"atoms.numbers" : 1, "output.total_energy" : 1, 
            "output.is_converged" : 1, "source_id" : 1, "operations" : 1,
            "nanoclusters" : 1
            }
        )
    simulations_list = [document for document in cursor]
    simulation_ids = [str(document["_id"]) for document in simulations_list]
    simulations  = dict(zip(simulation_ids, simulations_list))
    print("get simulations from elements")
    return simulations, simulations_list

def onpick3(event):
    ind = event.ind
    print('onpick3 scatter:', ind)
    try:
        first_id = ind[0]
    except:
        return
    #print(first_id)
    simulation_id = dft_ids[first_id]
    simulation = simulations[str(simulation_id)]
    pp(simulation)
    print("simulation id", simulation_id)
    atoms = atoms_dict_to_ase(simulation["atoms"])
    view(atoms)


def get_ml(workflow_id, is_test):
    if is_test:
        db_name = "testdb"
    else:
        db_name = "ncdb"

    ext_db, _ = get_external_database(
        username = "myusername",
        password = "mypassword",
        db_name = db_name,
        )
    ml = ext_db["machine_learning"].find({"workflow_id": workflow_id}, {"metrics_validation" : 1,"metrics_training" : 1,"metrics_test" : 1,"workflow_id" : 1,
        "training_set": 1, "test_set" : 1, "output" : 1, "method_params" : 1,
        })
    return ml


def get_property(workflow_id, calc_ids, is_test = False):
    if is_test:
        db_name = "testdb"
    else:
        db_name = "ncdb"

    ext_db, extdb_connect = get_external_database(
        username = "myusername",
        password = "mypassword",
        db_name = db_name,
        )

    simulations = fetch_simulations(extdb_connect, calc_ids)
    reaction_energies_list = np.zeros(len(calc_ids)).tolist()
    for idx, calc_id in enumerate(calc_ids):
        simulation = simulations[str(calc_id)]
        # get current simulation total_energy
        simulation_total_energy = simulation["output"].get("total_energy", 0.0)
        # iterate over
        # adsorbates
        adsorbates = simulation["adsorbates"]
        # nanoclusters
        nanoclusters = simulation["nanoclusters"]
        # substrates
        substrates = simulation["substrates"]

        component_types = [adsorbates, nanoclusters, substrates]
        reaction_energy = simulation_total_energy
        print("energy before adding references", reaction_energy)
        for components in component_types:
            for component in components:
                reference_id = component["reference_id"]
                print(reference_id)
                try:
                    reference_simulation = simulations[str(reference_id)]
                except:
                    print("getting reference from database")
                    reference_simulation  = ext_db["simulations"].find_one({"_id": reference_id})
                    simulations[str(reference_id)] = reference_simulation
                try:
                    total_energy = reference_simulation["output"]["total_energy"]
                except:
                    print("total_energy not found! Not contributing to reaction energy!")
                    total_energy = 0.0
                try:
                    reaction_energy -= float(total_energy)
                except:
                    print("Energy not understood!")
                    print(total_energy)

                print(reaction_energy, "reference", reference_id)
        reaction_energies_list[idx] = reaction_energy
    return reaction_energies_list, simulations


def get_learning_curve(workflow_id, is_test):
    # train, validation and test MAE
    training_sizes = []
    mae_train = []
    mae_validation = []
    mae_test = []

    if True:
        ml = get_ml(workflow_id, is_test)
        for entry in ml:
            train_size = len(entry["training_set"])
            training_sizes.append(train_size)
            test_size = len(entry["test_set"])
            print("train size", train_size, "test_size", test_size)
            mae_train.append(entry["metrics_training"]["mae"] * 27.2114)
            mae_validation.append(entry["metrics_validation"]["mae"] * 27.2114)
            mae_test.append(entry["metrics_test"]["mae"] * 27.2114)



    def draw_func():
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        ax.plot(training_sizes, mae_train, label = 'training set')
        ax.plot(training_sizes, mae_validation, label = 'validation set')
        ax.plot(training_sizes, mae_test, label = 'test set')
        ax.xaxis.label.set_fontsize(20)
        ax.set_xlabel('training set size', size=20)
        #ax.set_ylabel(r'$\alpha^2$ F [Ha]', size=20)
        ax.set_ylabel('MAE [eV]', size=20)
        ax.ticklabel_format(style='sci', axis='y', size =20, scilimits=(0,0))
        ax.legend(loc = 'upper right')
        plt.show()
        return None

    draw_func()
    return


def get_parity_plot(ml_document):
    global dft_ids, simulations

    # adsorption energies
    test_ids = ml_document["test_set"]
    training_ids = ml_document["training_set"]
    #prediction_ids = entry["prediction_set"]
    dft_ids = []
    dft_ids.extend(training_ids)
    dft_ids.extend(test_ids)

    test_property = ml_document["output"]["label_test"]
    training_property = ml_document["output"]["label_train"]
    predicted_property = ml_document["output"]["label_predicted"]
    predicted_property = []
    predicted_property.extend(training_property)
    predicted_property.extend(test_property)
    #y_property.extend(predicted_property)
    predicted_property = np.array(predicted_property)

    dft_property, simulations = get_property(ml_document["workflow_id"], dft_ids, is_test = False)
    dft_property = np.array(dft_property)


    # get_simulations_from_elements
    # overlap of ids
    # plot certain points in different colour

    target_element = np.zeros(len(dft_ids), dtype=bool)
    for idx, simulation in simulations.items():
        if simulation["atoms"]["numbers"].count(78) == 55:
            print("pure platinum")
            print(idx)
            target_element[int(idx) == np.array(dft_ids).astype(int)] = 1

    print("pure platinum", target_element.sum())


    matplotlib.rcParams.update({'font.size': 18})

    fig, (ax0) = plt.subplots(1, 1)
    fig.canvas.mpl_connect('pick_event', onpick3)

    # convert to eV
    x = dft_property * 27.2114
    y = predicted_property * 27.2114

    print(x.shape, y.shape)

    plt.scatter(x,y, picker = True)
    plt.scatter(x[target_element],y[target_element], c = 'r')


    x_line = np.linspace(ax0.get_xlim(), ax0.get_ylim(), 100)
    lims = [
        np.min([ax0.get_xlim(), ax0.get_ylim()]),  # min of both axes
        np.max([ax0.get_xlim(), ax0.get_ylim()]),  # max of both axes
    ]
    # now plot both limits against eachother
    #ax0.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax0.plot(x_line, x_line, '-r', zorder = 0)
    ax0.plot(x_line, x_line + 0.1, ':b')
    ax0.plot(x_line, x_line - 0.1, ':b')
    ax0.set_aspect('equal')
    ax0.set_xlim(lims)
    ax0.set_ylim(lims)
    ax0.set_xlabel(r'Calculated $\Delta$E [eV]', size=20)
    #ax0.set_ylabel(r'$\alpha^2$ F [Ha]', size=20)
    ax0.set_ylabel(r'Predicted $\Delta$E [eV]', size=20)
    plt.show()


    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('ids', metavar='N', type=int, nargs='+',
                    help='workflow id')
    parser.add_argument('--test', action='store_true')

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
    get_parity_plot(entry)



    get_learning_curve(args.ids[0], is_test)


