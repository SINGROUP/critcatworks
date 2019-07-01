# functions to store data in an external MongoDB database
import pymongo
from pprint import pprint as pp

def get_external_database(extdb_connect):
    """A helper function to connect to a mongodb database.

    Args:
        extdb_connect (dict):   dictionary containing the keys host,
                                username, password, authsource and db_name.

    Returns:
        pymongo object : address to database
    """
    CLIENT = pymongo.MongoClient(extdb_connect["host"], 
        username = extdb_connect["username"],
        password = extdb_connect["password"],
        authSource = extdb_connect["authsource"])
    db = CLIENT[extdb_connect["db_name"]]
    return db


def _reset_IDs_collection(db = None):
    """Helper function to reset the counting of the following ids:
    simulations
    workflows
    machine_learning

    Args:
        db (pymongo object) : address to database

    Returns:
        int : number of documents in the IDs collection (expectedly 1)
    """
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


def _query_id_counter_and_increment(collection, db):
    """Helper function to query an ID and increase it by one.

    Args:
        collection (str) :  name of the collection id to fetch and increment.
                            Valid strings are 'simulations', 'workflows' and
                            'machine_learning'
        db (pymongo object) : address to database

    Returns:
        int : id counter of the specified collection.
    """
    ids_collection = db['IDs']
    id_counter = ids_collection.find_one_and_update({}, {'$inc': {collection: 1}})[collection]
    return id_counter

def update_simulations_collection(extdb_connect, **kwargs):
    """A new document is added to the simulations collection of the mongodb database.
    It contains records of all manipulation steps of a structure, in particular the 
    initial structure, structure after DFT relaxation, 
    structure with added or removed asdorbates, etc.
    The documents should be in a specific format. Any arguments can be
    specified, however, the optional arguments below should be 
    consistently given to allow for comprehensive database querying.
                            
    Args:
        extdb_connect (dict):   dictionary containing the keys host,
                                username, password, authsource and db_name. 

        source_id (int)     : ID of the parent simulation that originated this, -1 if none
        workflow_id (int)   : ID of workflow when instance was added, -1 if none
        wf_sim_id (int)     : ID of simulation (unique within the workflow this belongs to)
        atoms (dict)        :   dictionary with information about the atoms.
                                should be in the following format

                                numbers (1D ndarray)   : list of atomic numbers as numpy array [N] of ints
                                positions (2D ndarray) : positions as numpy matrix [Nx3] of doubles
                                constraints (2D ndarray) : frozen flags a matrix [Nx3] of int [optional] 1 = frozen, 0 = free
                                pbc (bool)             : use periodic boundaries
                                cell (2D ndarray)      : matrix 3x3 with cell vectors on the rows
                                celldisp (1D ndarray)  : displacement of cell from origin
                                info (dict)            : field for additional information related to structure
        nanoclusters (list of ATOMS dict) :   list of dictionaries with information about the nanocluster(s)
                                        The dictionaries should have the following form:

                                        reference_id (int) : ID of the simulation where this cluster was made, -1 if original
                                        atom_ids (1D ndarray) : atom indices in the ATOMS dictionary of the simulation record.

        adsorbates (list of dict) :     list of dictionaries with information about the adsorbate(s)
                                        The dictionaries should have the following form:

                                        reference_id (int) : ID of the simulation to use as reference
                                        atom_ids (1D ndarray) : atom indices in the ATOMS dictionary of the simulation record.
                                        site_class (str) : class of adsorption site: “top”, “bridge”, “hollow”, “4-fold hollow”
                                        site_ids (1D ndarray) : list of atom ids (in simulation record) that define the adsorption site

        substrate (list of dict) :      list of dictionaries with information about the substrate(s)
                                        The dictionaries should have the following form:

                                        reference_id (int) : ID of the parent support simulation, -1 if no parent
                                        atom_ids (1D ndarray) : atom indices in the corresponding ATOMS dictionary

        operations (list) : List of dictionaries, each describing one operation. Always with respect to the parent simulation if applicable.
                            The dictionaries can be of arbitrary form.
        inp (dict)        : property/value pairs describing the simulation input
                            The dictionary can be of arbitrary form.
        output (dict)     : property/value pairs output by the calculation
                            The dictionary can be of arbitrary form.

    Returns:
        dict :  A dictionary with the provided arguments plus
                a unique id provided by the database.
    """
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
    """A new document is added to the workflows collection of the mongodb database.
    (Usually at the beginning of the workflow run.)
    It contains records of all types of workflows. The documents should be in a specific format. 
    Any arguments can be specified, however, certain arguments below should be 
    consistently given to allow for comprehensive database querying.
                            
    Args:
        extdb_connect (dict):   dictionary containing the keys host,
                                username, password, authsource and db_name. 
        username (str) :        user who executed the workflow
        creation_time (str) :   time of creation of the workflow
        parameters (dict) :     workflow-specific parameters
        name (str) :            custom name of workflow
        workflow_type (str) :   custom type of workflow

    Returns:
        dict :  Contains the keys username, name, workflow_type, creation_time,
                parameters and _id, the latter being a 
                unique id provided by the database.
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
    """A new document is added to the machine_learning collection of the mongodb database.
    It contains records of all types of workflows. The documents should be in a specific format. 
    Any arguments can be specified, however, certain arguments below should be 
    consistently given to allow for comprehensive database querying.

    Args:
        workflow_id (int) : ID of workflow which the machine learning run was part of
        method (str) : name of the ML method: krr, nn, ...
        method_params (dict) : Parameters of the method
        descriptor (str) : name of the descriptor: soap, mbtr, lmbtr, cm, ...
        descriptor_params (dict) : Parameters of the descriptor used
        training_set (1D ndarray) : list of simulation IDs used for training
        validation_set (1D ndarray) :   list of simulation IDs used in validation. 
                                        If empty, cross-validation was used.
        test_set (1D ndarray) : list of simulation IDs used in testing. 
                                If empty, only validation was used
        prediction_set (1D ndarray) : list of simulation IDs used for prediction.
        metrics_training (dict) : dictionary of (“metric name”: value) on training set
                                key (str) : name of the metric
                                value (float) : calculated value
        metrics_validation (dict) : dictionary of (“metric name”: value) on validation set
                                    key (str) : name of the metric
                                    value (float) : calculated value

        metrics_test (dict) :   dictionary of (“metric name”: value) on test set
                                key (str) : name of the metric
                                value (float) : calculated value
        output (dict)       :   relevant training output info

    Returns:
        dict :  A dictionary with the provided arguments plus
                a unique id provided by the database.
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
    """Fetches simulation records by simulation id

    Args:
        extdb_connect (dict):   dictionary containing the keys host,
                                username, password, authsource and db_name. 
        simulation_ids (1D ndarray) : unique identifiers of the simulation collection.

    Returns:
        list : documents of the simulation collection
    """
    db = get_external_database(extdb_connect)
    cursor =   db['simulations'].find({"_id" : {"$in" : simulation_ids }})
    simulations_list = [document for document in cursor]
    simulation_ids = [str(document["_id"]) for document in simulations_list]
    simulations  = dict(zip(simulation_ids, simulations_list))
    print(simulation_ids, len(simulation_ids))
    return simulations


def gather_all_atom_types(calc_ids, simulations):
    """Helper function to determine all atom types in the dataset

    Args:
        calc_ids (list) : ids of the simulation collection
        simulations (list) : simulation documents

    Returns:
        list :  a sorted unique list of atomic numbers in the
                dataset
    """
    # going through nc atoms once to find atom types
    atomic_numbers = []
    for idx, calc_id in enumerate(calc_ids):
        atoms_dict = simulations[str(calc_id)]["atoms"]
        atoms = atoms_dict_to_ase(atoms_dict)
        atomic_numbers.extend(atoms.get_atomic_numbers())

    sorted_list_atomic_numbers = list(sorted(set(atomic_numbers)))

    all_atomtypes = sorted_list_atomic_numbers
    return all_atomtypes
