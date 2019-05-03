# functions to store data in an external MongoDB database
import pymongo
from pprint import pprint as pp

def get_external_database(extdb_connect):
    CLIENT = pymongo.MongoClient(extdb_connect["host"], 
        username = extdb_connect["username"],
        password = extdb_connect["password"],
        authSource = extdb_connect["authsource"])
    db = CLIENT[extdb_connect["db_name"]]
    return db


# helper function to reset entry in collection IDs
def _reset_IDs_collection(db = None):

    # counting collection
    ids_collection = db['IDs']
    #print(id_collection)
    starting_ids = {'_id': -1,
            'simulations' : 1,
            'workflows' : 1,
            'machine_learning' : 1,
            }

    count = ids_collection.count({})
    if count == 0:
        result = ids_collection.insert_one(starting_ids)
        print("no entry yet")
        print(result)
    else:
        ids_collection.update({"_id" : -1}, starting_ids)
        print(list(db['IDs'].find()))
    return count


# helper function to query an ID and increase it by one
def _query_id_counter_and_increment(collection, db):
    ids_collection = db['IDs']
    id_counter = ids_collection.find_one_and_update({}, {'$inc': {collection: 1}})[collection]
    return id_counter


# update simulation
#def update_simulations_collection(atoms = {}, 
#    source_id = -1, workflow_id = -1, 
#    nanoclusters = [], adsorbates = [], substrates = [], 
#    operations = [], inp = {}, output = {},
#    **kwargs):
def update_simulations_collection(extdb_connect, **kwargs):
    """
    # generated cluster manual
    # and automated using clusterer
    # generated adsorption site
    # generated coverage (added/removed adsorbates)
    # DFT simulation
    """
    # construct dictionary
    #dct = {
    #    "atoms" : atoms, 
    #    "source_id" : source_id,
    #    "workflow_id" : workflow_id,
    #    "nanoclusters" : nanoclusters,
    #    "adsorbates" : adsorbates,
    #    "substrates" : substrates,
    #    "operations" : operations,
    #    "inp" : inp,
    #    "output" : output,
    #    }
    #dct.update(kwargs)

    dct = kwargs

    db = get_external_database(extdb_connect)
    simulations = db['simulations']
    # request id counter
    simulation_id = _query_id_counter_and_increment('simulations', db)

    dct['_id'] = simulation_id
    print("new simulation id: " , str(simulation_id))

    simulations.insert_one(dct)
    return dct

# workflows

def update_workflows_collection(username, password, creation_time, 
    extdb_connect, parameters = {},
    name = "UNNAMED", workflow_type = "NO_TYPE",
    **kwargs):
    """
    Update after workflow has finished
    """
    db = get_external_database(extdb_connect)
    workflows = db['workflows']

    # request id counter
    workflow_id = _query_id_counter_and_increment('workflows', db)

    # construct dictionary
    dct = {'_id' : workflow_id,
        "username" : username,
        "name" : name,
        "workflow_type" : workflow_type,
        "creation_time" : creation_time,
        "parameters" : parameters,
        }
    dct.update(kwargs)

    workflows.insert_one(dct)
    return dct

# machine_learning

def update_machine_learning_collection(method, extdb_connect, workflow_id = -1, 
    method_params = {}, descriptor = "soap",
    descriptor_params = {},
    training_set = [], validation_set = [], test_set = [],prediction_set = [],
    metrics_training = {}, metrics_validation = {}, metrics_test = {},
    output = {},
    **kwargs):
    """
    update after each machine learning run
    """
    db = get_external_database(extdb_connect)
    machine_learning = db['machine_learning']

    # request id counter
    machine_learning_id = _query_id_counter_and_increment('machine_learning', 
        db)

    # construct dictionary
    dct = {'_id' : machine_learning_id,
        "method" : method,
        "workflow_id" : workflow_id,
        "method_params" : method_params,
        "descriptor" : descriptor,
        "descriptor_params" : descriptor_params,
        "training_set" : training_set,
        "validation_set" : validation_set,
        "test_set" : test_set,
        "prediction_set" : prediction_set,
        "metrics_training" : metrics_training,
        "metrics_validation" : metrics_validation,
        "metrics_test" : metrics_test,
        "output" : output,
        }
    dct.update(kwargs)

    machine_learning.insert_one(dct)
    return dct

def fetch_simulations(extdb_connect, simulation_ids):
    db = get_external_database(extdb_connect)
    cursor =   db['simulations'].find({"_id" : {"$in" : simulation_ids }})
    simulations_list = [document for document in cursor]
    simulation_ids = [str(document["_id"]) for document in simulations_list]
    simulations  = dict(zip(simulation_ids, simulations_list))
    print(simulation_ids, len(simulation_ids))
    return simulations


### custom types conversion functions

if __name__ == "__main__":
    #reset = _reset_IDs_collection()
    #print(reset)

    #id_counter = _query_id_counter_and_increment("simulations")
    #print(id_counter)
    ids = [2222,2500, 3333, 2225, 2000, 2400, 3001]
    simulations = fetch_simulations({"host" : "nanolayers.dyndns.org:27017", "username" : "mjcritcat", "password" : "heterogeniuscatalysis", "db_name" : "testdb", "authsource" : "testdb"}, ids)
    #print(len(simulations))
    #for document in cursor:
    #    print(document)
    print(len(simulations))
    print(simulations.keys())
    for idx in ids:
        print(simulations[str(idx)]["_id"], "==", idx)


    #wfc = update_workflows_collection("random_dude", "Sunday")

    #sc = update_simulations_collection(atoms = "HERE_SHOULD_BE_A_CUSTOM_TYPE", notes = "adding_keys_allowed")

    #mlc = update_machine_learning_collection("krr", 
    #    russian_exchange_student = "Sakmiov")

    #print(sc)
    #print(wfc)
    #print(mlc)
