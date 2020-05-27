# tool to visualize a structure from a mongodb database 
# requires the id(s) to be visualized
# default is the critcat database, the  flag --is_test redirects to the test database
# usage: python3 visualize_structure 250 251
# 

#from critcatworks.database.format import atoms_dict_to_ase, ase_to_atoms_dict
# removed critcatworks dependency


def atoms_dict_to_ase(atoms_dict):
    """
    Helper function to convert a ATOMS dictionary into an 
    ase.Atoms object
    Args:
        atoms_dict (dict) : dictionary with information about the atoms.
                            should be in the following format
                            numbers (1D ndarray)   : list of atomic numbers as numpy array [N] of ints
                            positions (2D ndarray) : positions as numpy matrix [Nx3] of doubles
                            constraints (2D ndarray) : frozen flags a matrix [Nx3] of int [optional] 1 = frozen, 0 = free
                            pbc (bool)             : use periodic boundaries
                            cell (2D ndarray)      : matrix 3x3 with cell vectors on the rows
                            celldisp (1D ndarray)  : displacement of cell from origin
                            info (dict)            : field for additional information related to structure
    Returns:
        ase.Atoms : Corresponding ase.Atoms object
    """
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


def ase_to_atoms_dict(atoms):
    """
    Helper function to convert an 
    ase.Atoms object into its corresponding python dictionary
    Args:
        atoms (ase.Atoms) : ase.Atoms object
    Returns:
        dict : Corresponding python dictionary
    """
    positions = atoms.get_positions().tolist()
    cell = atoms.get_cell().tolist()
    pbc = atoms.get_pbc().tolist()
    numbers = atoms.get_atomic_numbers().tolist()
    try:
        constraints = atoms.constraints.tolist()
    except:
        constraints = atoms.constraints
    celldisp = atoms.get_celldisp().tolist()
    info = atoms.info

    atoms_dict = {
        "positions" : positions,
        "cell" : cell,
        "pbc" : pbc,
        "numbers" : numbers,
        "constraints" : constraints,
        "celldisp" : celldisp,
        "info" : info,    
        }
    return atoms_dict


import ase, ase.io
from ase.visualize import view
import pymongo
from pprint import pprint as pp
import argparse
import numpy as np


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



def get_simulation(simulation_id, is_test):
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('ids', metavar='N', type=int, nargs='+',
                    help='simulation id')
    parser.add_argument('--test', action='store_true')

    args = parser.parse_args()
    print(args.ids)
    is_test = args.test
    for idx in args.ids:
        simulation = get_simulation(idx, is_test)
        pp(simulation)
        print("positions:", np.array(simulation["atoms"]["positions"]).shape)
        atoms = atoms_dict_to_ase(simulation["atoms"])
        ase.io.write("testout_" + str(idx) + ".xyz", atoms)
        view(atoms)
