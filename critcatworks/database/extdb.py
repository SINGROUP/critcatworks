# functions to store data in an external MongoDB database
import pymongo

def get_external_database(host = "mongodb+srv://austerity-hgeov.mongodb.net/test", 
    database = "test",username = "marc", password = 'marcrulez0r'):
    CLIENT = pymongo.MongoClient(host, username = username, password = password)
    db = CLIENT[database]
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
def update_simulations_collection(wf_sim_id, 
    atoms, 
    source_id = -1, workflow_id = -1, 
    nanoclusters = [], adsorbates = [], substrates = [], 
    operations = [], inp = {}, output = {},
    **kwargs):
    """
    IS wf_sim_id really necessary? update ids in internal database 
    after querying?
    
    # generated cluster manual
    # and automated using clusterer
    # generated adsorption site
    # generated coverage (added/removed adsorbates)
    # DFT simulation
    """
    db = get_external_database()
    simulations = db['simulations']
    # request id counter
    simulation_id = _query_id_counter_and_increment('simulations', db)



    # construct dictionary
    dct = {'_id' : simulation_id,
        "wf_sim_id" : wf_sim_id,
        "atoms" : atoms, 
        "source_id" : source_id,
        "workflow_id" : workflow_id,
        "nanoclusters" : nanoclusters,
        "adsorbates" : adsorbates,
        "substrates" : substrates,
        "operations" : operations,
        "inp" : inp,
        "output" : output,
        }
    dct.update(kwargs)

    simulations.insert_one(dct)
    return dct

# workflows

def update_workflows_collection(username, creation_time, 
    parameters = {},
    name = "UNNAMED", workflow_type = "NO_TYPE",
    **kwargs):
    """
    Update after workflow has finished
    """
    db = get_external_database()
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

def update_machine_learning_collection(method, workflow_id = -1, 
    method_params = {}, descriptor = "soap",
    descriptor_params = {},
    training_set = [], validation_set = [],
    metrics_training = {}, metrics_validation = {},
    output = {},
    **kwargs):
    """
    update after each machine learning run
    """
    db = get_external_database()
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
        "metrics_training" : metrics_training,
        "metrics_validation" : metrics_validation,
        "output" : output,
        }
    dct.update(kwargs)

    machine_learning.insert_one(dct)
    return dct

### custom types conversion functions

if __name__ == "__main__":
    #reset = _reset_IDs_collection()
    #print(reset)

    id_counter = _query_id_counter_and_increment("simulations")
    print(id_counter)

    wfc = update_workflows_collection("random_dude", "Sunday")

    sc = update_simulations_collection(-42, 
        "HERE_SHOULD_BE_A_CUSTOM_TYPE", notes = "adding_keys_allowed")

    mlc = update_machine_learning_collection("krr", 
        russian_exchange_student = "Sakmiov")

    print(sc)
    print(wfc)
    print(mlc)
