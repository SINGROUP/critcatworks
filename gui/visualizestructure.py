from critcatworks.database.format import atoms_dict_to_ase, ase_to_atoms_dict
import ase, ase.io
from ase.visualize import view
import pymongo

def get_external_database(host = "mongodb+srv://austerity-hgeov.mongodb.net/test", 
    database = "test",username = "marc", password = 'marcrulez0r'):
    CLIENT = pymongo.MongoClient(host, username = username, password = password)
    db = CLIENT[database]
    return db

def get_simulation(simulation_id):
    ext_db = get_external_database()
    simulation = ext_db["simulations"].find_one({"_id": simulation_id})
    return simulation

if __name__ == '__main__':
    ID = 2149
    simulation = get_simulation(ID)

    atoms = atoms_dict_to_ase(simulation["atoms"])

    view(atoms)
