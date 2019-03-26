from critcatworks.database.format import atoms_dict_to_ase, ase_to_atoms_dict
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



def get_simulation(simulation_id):
    ext_db = get_external_database(
        username = "mjcritcat",
        password = "heterogeniuscatalysis",
        #db_name = "testdb",
        db_name = "ncdb",
        )
    simulation = ext_db["simulations"].find_one({"_id": simulation_id})
    return simulation

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('ids', metavar='N', type=int, nargs='+',
                    help='simulation id')

    args = parser.parse_args()
    print(args.ids)
    for idx in args.ids:
        simulation = get_simulation(idx)
        pp(simulation)
        print("positions:", np.array(simulation["atoms"]["positions"]).shape)
        atoms = atoms_dict_to_ase(simulation["atoms"])
        view(atoms)
