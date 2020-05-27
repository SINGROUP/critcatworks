# utility function for the nanocluster workflow
# ATTENTION: Some parameters have to be modified in this script,
# such as: compositions, elements, natoms (nanocluster size)
# given a workflow_id
# plots the excess energies (w.r.t pure element nanoclusters)

import pymongo, argparse
import numpy as np
import ase
from matplotlib import pyplot as plt
from collections import defaultdict

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
    return db



def get_data_from_wflow(workflow_id = 45, is_test = False):
    if is_test:
        db_name = "testdb"
    else:
        db_name = "ncdb"

    ext_db = get_external_database(
        username = "myusername",
        password = "mypassword",
        db_name = db_name,
        )
    data = ext_db.simulations.find(
        {'workflow_id' : workflow_id, 
        }, 
        {
            "atoms.numbers" : 1, 
            "output.total_energy" : 1, 
            "output.is_converged" : 1, 
            "source_id" : 1,
            "inp.cEA" : 1,
            "inp.cEB" : 1,
            "inp.eAA" : 1,
            "inp.eAB" : 1,
            "inp.eBB" : 1,
            }
        )
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('ids', metavar='N', type=int, nargs='+',
                    help='workflow id')
    parser.add_argument('--test', action='store_true')

    # input:    wf id
    args = parser.parse_args()
    print(args.ids)
    is_test = args.test
    workflow_id = args.ids[0]

    natoms = 55
    compositions = [6, 13, 28, 42, 49]
    elements = ["Fe","Co", "Ni","Cu","Ti", "Pt"]

    # load dataset from nanocluster workflow
    database = get_data_from_wflow(is_test = False)
    database = list(database)

    # prepare subplots
    fig, axes = plt.subplots(len(elements) -1 , len(elements) -1, sharex='col', sharey = 'row',
        gridspec_kw={'wspace': 0, 'hspace' : 0.1})
    #fig.set_size_inches(11.69,8.27)
    fig.set_size_inches(11.69 * 1.5, 16.53 * 1.5)
    plt.setp(axes, xticks=[compositions[0], compositions[2], compositions[4]], 
        #xticklabels=['a', 'b', 'c'],
        #yticks=[1, 2, 3]
        )

    # element combinations, compositions and 
    # entries per composition
    for el1_idx, element1 in enumerate(elements):
        atomic_number_1 = ase.data.atomic_numbers[element1]
        for el2_idx, element2 in enumerate(elements):
            atomic_number_2 = ase.data.atomic_numbers[element2]
            if el1_idx >= el2_idx:
                pass
            else:
                simulations = [x for x in database if set(x['atoms']['numbers']).issubset({atomic_number_1, atomic_number_2})]
                # get pure nanoclusters
                pure_sim1_lst = [x for x in database if set(x['atoms']['numbers']).issubset({atomic_number_1})]
                pure_sim2_lst = [x for x in database if set(x['atoms']['numbers']).issubset({atomic_number_2})]
                print("amount of pure simulations", len(pure_sim1_lst), len(pure_sim2_lst))
                for pure_sim1, in zip(pure_sim1_lst, ):
                    try:
                        total_energy = pure_sim1_lst["output"]["total_energy"]
                        break
                    except:
                        continue
                for pure_sim2 in pure_sim2_lst:
                    try:
                        total_energy = pure_sim2_lst["output"]["total_energy"]
                        break
                    except:
                        continue

                
                #print(pure_sim1, pure_sim2)

                x_compositions = []
                ee_lst = []

                print("initial query of", len(simulations), "documents")
                for comp_idx, composition in enumerate(compositions):
                    simulation_ids = []
                    source_ids = []
                    en_dict = {}
                    core_dict = {}
                    segregation_dict = {}
                    c_value_dict = {}

                    for simulation in simulations:
                        if (int(simulation["atoms"]["numbers"].count(atomic_number_1)) != composition):
                            # sort by composition
                            pass
                        else:
                            simulation_ids.append(simulation["_id"])
                            # get energies
                            try:
                                is_converged = simulation["output"]["is_converged"]
                            except:
                                is_converged = False
                            if not is_converged:
                                #print("not converged")
                                pass
                            else:
                                # get total energies
                                total_energy = simulation["output"]["total_energy"]

                                # compute excess energy
                                # [E(mix = n * el1 + (m-n) * el2) - n/m * E(pure1) - (m-n)/m E(pure2) ]/ m
                                # m = number of atoms, n = number of atoms of pure element 1

                                pure_toten1 = pure_sim1["output"]["total_energy"]
                                pure_toten2 = pure_sim2["output"]["total_energy"]
                                ee = (total_energy - composition / natoms * pure_toten1 - (natoms - composition) / natoms * pure_toten2 ) / natoms
                                ee = ee * 27.2114
                                #print("excess energy:", ee)
                                x_compositions.append(composition)
                                ee_lst.append(ee)
                # add subplot
                ax = axes[el1_idx, el2_idx -1]
                #x_ratio = np.array(x_compositions) / natoms
                ax.scatter(x_compositions, ee_lst, marker='o', facecolors='none', edgecolors = 'b')
                ax.set_title(r'${}_n {}_{{{}}}$'.format(element1, element2, '55-n'))
                

                comp_ee_dct = defaultdict(lambda:0)
                for ee_i, ee in enumerate(ee_lst):
                    if comp_ee_dct[x_compositions[ee_i]] < ee:
                        pass
                    else:
                        comp_ee_dct[x_compositions[ee_i]] = ee

                xy_pareto = sorted(comp_ee_dct.items())
                x_pareto = [c for c,_ in xy_pareto]
                y_pareto = [e for _,e in xy_pareto]
                ax.plot(x_pareto, y_pareto, '--r')
                ax.axhline(ls ='--', color = 'black')

# plot excess energies
fig.text(0.5, 0.08, 'n', ha='center')
fig.text(0.04, 0.5, 'excess energy [eV / atom]', va='center', rotation='vertical')

fig.savefig('nc_excess_energies.png', bbox_inches = 'tight')
#plt.show()
