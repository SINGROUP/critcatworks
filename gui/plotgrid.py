# utility script for nanoclusters workflow
# offline or database
# plotting in 2D color grid the energies of 
# each binary TM structure
# interactive, view structures upon clicking

import argparse
import pymongo
import ase, ase.io
from ase.visualize import view
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.pyplot import colorbar
import time

def atoms_dict_to_ase(atoms_dict):
    cell = atoms_dict['cell']
    celldisp = atoms_dict['celldisp']
    constraints = atoms_dict['constraints']
    pbc = atoms_dict['pbc']
    numbers = atoms_dict['numbers']
    positions = atoms_dict['positions']
    info = atoms_dict.get("info", {})

    atoms = ase.Atoms(numbers=numbers,
        positions=positions,
        cell = cell,
        celldisp = celldisp,
        constraint = constraints,
        pbc = pbc,
        info = info)
    return atoms

def onpick3(event):
    ind = event.ind
    print('onpick3 scatter:', ind, 
        np.take(Z, ind), 
        np.take(Z, ind))
    first_id = ind[0]
    print(first_id)

def onclick(event):
    # x and y seem to be exchanged
    y, x = event.xdata, event.ydata
    int_x, int_y = int(x), int(y)

    simulation_id = sim_id_matrix[int_x, int_y]
    if simulation_id == 0:
        pass
    else:
        simulation = get_simulation(simulation_id)
        atoms = atoms_dict_to_ase(simulation["atoms"])
        view(atoms)
    print(simulation_id)


class Cursor(object):
    def __init__(self, ax):
        self.ax = ax
        self.lx = ax.axhline(color='k')  # the horiz line
        self.ly = ax.axvline(color='k')  # the vert line

        # text location in axes coords
        self.txt = ax.text(0.7, 0.9, '', transform=ax.transAxes)
        self.formula = None

    def mouse_move(self, event):
        if not event.inaxes:
            return

        x, y = event.xdata, event.ydata
        # x and y seem to be exchanged
        int_x, int_y = int(y), int(x)
        # update the line positions
        self.lx.set_ydata(y)
        self.ly.set_xdata(x)
        try:
            formula = self.formula[int_x, int_y]
        except:
            formula = ""

        self.txt.set_text(formula)
        self.ax.figure.canvas.draw()

# access to database
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


def get_simulation(simulation_id, is_test = False):
    if is_test:
        db_name = "testdb"
    else:
        db_name = "ncdb"

    ext_db = get_external_database(
        username = "myusername",
        password = "mypassword",
        db_name = db_name,
        )
    simulation = ext_db["simulations"].find_one({"_id": simulation_id})
    return simulation


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
    simulations = ext_db.simulations.find(
        {'workflow_id' : workflow_id, 
        "$expr":{"$setIsSubset":["$atoms.numbers",[el1,el2]]}
        }, 
        {"atoms.numbers" : 1, "output.total_energy" : 1, "output.is_converged" : 1, "source_id" : 1}
        )

    return simulations


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


a = time.time()
# offline mode
# argparse
parser = argparse.ArgumentParser(description='Process some integers.')
#parser.add_argument('ids', metavar='N', type=int, nargs='+',
#                help='simulation id')
parser.add_argument('--offline', action='store_true')
parser.add_argument('--core', action='store_true')
parser.add_argument('--segregation', action='store_true')

args = parser.parse_args()
is_offline = args.offline
is_core = args.core
is_segregation = args.segregation


# loop over structures
natoms = 55
n_entries = 10
compositions = [6, 13, 28, 42, 49]
elements = ["Fe","Co", "Ni","Cu","Ti", "Pt"]
x_dimensions = len(elements) * len(compositions)
y_dimensions = len(elements) * n_entries
Z = np.zeros((x_dimensions, y_dimensions))
print("property grid shape", Z.shape)
formula = np.empty([x_dimensions, y_dimensions], dtype="U10")
sim_id_matrix = np.zeros((x_dimensions, y_dimensions))
#print(formula.shape)

if not is_offline:
    database = get_data_from_wflow(is_test = False)
    database = list(database)

    # element combinations, compositions and 
    # entries per composition
    for el1_idx, element1 in enumerate(elements):
        atomic_number_1 = ase.data.atomic_numbers[element1]
        for el2_idx, element2 in enumerate(elements):
            atomic_number_2 = ase.data.atomic_numbers[element2]
            if el1_idx > el2_idx:
                pass
            else:
                simulations = [x for x in database if set(x['atoms']['numbers']).issubset({atomic_number_1, atomic_number_2})]
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
                            pass
                            #print(int(simulation["atoms"]["numbers"].count(atomic_number_1)))
                        else:
                            try:
                                source_ids.append(simulation['source_id'])
                            except:
                                pass
                            simulation_ids.append(simulation["_id"])
                            # get energies (atomization energies?, cohesive energies?)
                            try:
                                is_converged = simulation["output"]["is_converged"]
                                if not is_converged:
                                    c_value = 0.3
                                else:
                                    total_energy = simulation["output"]["total_energy"]
                                    en_dict[simulation["_id"]] = total_energy
                                    # coreness metric
                                    cEA = simulation["inp"]["cEA"]
                                    cEB = simulation["inp"]["cEB"]
                                    core_dict[simulation["_id"]] = float(cEA) - float(cEB)
                                    # segregation metric
                                    eAA = simulation["inp"]["eAA"]
                                    eAB = simulation["inp"]["eAB"]
                                    eBB = simulation["inp"]["eBB"]
                                    segregation_dict[simulation["_id"]] = float(eAB) - 0.5 * (float(eAA) + float(eBB))
                                    c_value = 1.0
                            except:
                                c_value = 0.3
                            c_value_dict[simulation["_id"]] = c_value


                    simulation_ids = np.setdiff1d(simulation_ids, source_ids, assume_unique=True)
                    c_value_dict = {k:c_value_dict[k] for k in simulation_ids if k in c_value_dict}

                    if en_dict == {}:
                        pass
                    else:
                        x = en_dict.values()
                        xmin, xmax = min(x), max(x)
                        for key, value in en_dict.items():
                            if xmax == xmin:
                                c_value_dict[key] = 0.5
                            else:
                                c_value_dict[key] = (value-xmin)/(xmax-xmin) + 0.5
                    
                    # sort by energy
                    sorted_simulation_ids = sorted(c_value_dict.items(), key=lambda kv: kv[1])
                    sorted_idx = [key for key, value in sorted_simulation_ids]
                    simulation_ids = sorted_idx

                    # replace by other property
                    if is_core:
                        property_dct = core_dict
                    elif is_segregation:
                        property_dct = segregation_dict
                    else:
                        pass
                    if (is_core or is_segregation):
                        if property_dct == {}:
                            pass
                        else:
                            x = property_dct.values()
                            #xmin, xmax = min(x), max(x)
                            for key, value in property_dct.items():
                                c_value_dict[key] = value + 0.5
                                #c_value_dict[key] = (value-xmin)/(xmax-xmin) + 0.5

                    for entry in range(n_entries):
                        try:
                            simulation_id = simulation_ids[entry]
                            c_value = c_value_dict[simulation_id]
                        except:
                            c_value = 0
                            simulation_id = 0

                        x_id = el1_idx * len(compositions) + comp_idx
                        y_id = el2_idx * n_entries + entry

                        Z[x_id, y_id] = c_value
                        formula[x_id, y_id] = element1 + str(composition) + element2 + str(natoms - composition)
                        sim_id_matrix[x_id, y_id] = simulation_id
                    #print("final", len(simulation_ids))
    # save numpy arrays for offline use
    np.save("Z.npy", Z)
    np.save("formula.npy", formula)
    np.save("simulation_ids.npy", sim_id_matrix)
else:
    Z = np.load("Z.npy")
    formula = np.load("formula.npy")
    sim_id_matrix = np.load("simulation_ids.npy")

# plotting
b = time.time()
print("time", b -a)


fig, (ax0) = plt.subplots(1, 1)
c = ax0.pcolor(Z, picker=True)

ax0.set_xticklabels(elements + [""] )
ax0.set_yticklabels(elements + [""] )

print(compositions * len(elements))

cursor = Cursor(ax0)
cursor.formula = formula
fig.canvas.mpl_connect('motion_notify_event', cursor.mouse_move)
fig.canvas.mpl_connect('pick_event', onpick3)
cid = fig.canvas.mpl_connect('button_press_event', onclick)

plt.colorbar(c,ax=ax0)

#fig.tight_layout()

if is_core:
    prefix = "core_"
elif is_segregation:
    prefix = "segregation_"
else:
    prefix = "energy_"

plt.savefig(prefix + "test.png")
plt.show()



# draw a new figure and replot the colorbar there
fig,ax = plt.subplots(figsize=(2,3))
plt.colorbar(c,ax=ax)
ax.remove()
#plt.savefig('plot_onlycbar.png')
# save the same figure with some approximate autocropping
plt.savefig(prefix + 'plot_onlycbar_tight.png',bbox_inches='tight')




if __name__ == '__main__':
    pass

