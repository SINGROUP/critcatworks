.. _database:

Database
========

For storing information about the workflows, machine learning, 
DFT simlations and other structure manipulations, MongoDB 
is used as a permament, shared database.
MongoDB is a NoSQL database, with bson dictionaries as documents,
which is very close to python documents. If you know python,
MongoDB documents are easy to read and query.

Database format
---------------

In the following, we list the entries that data records need to have in each collection. A data record is assumed to be a python dictionary, and the tables give the keys that are expected to be found therein. The following collections are defined: IDs, simulations, workflows and machine_learning.


Collection: IDs

contains the ID values that will be assigned to the newest entries in other collections


_id (int):
    automatic compulsory internal id for MongoDB

simulations (int):
    ID to be used by the next new entry of type simulation

workflows (int):
    ID to be used by the next new entry of type workflows 

machine_learning (int):
    ID to be used by the next new entry of type machine_learning 

(notice how the field name matches the collection name)



Collection: simulations


_id (int): 
    unique identifier

source_id (int):
    ID of the parent simulation that originated this, -1 if none

workflow_id (int):
    ID of workflow when instance was added, -1 if none

wf_sim_id (int):
    ID of simulation (unique within the workflow this belongs to)

atoms (ATOMS):
    dictionary with information about the atoms.

nanoclusters (NANOCLUSTER):
    list of dictionaries with information about the nanocluster(s)
adsorbates (ADSORBATE):
    list of dictionaries with information about the adsorbate(s)

substrates (SUBSTRATE):
    list of dictionaries with information about the substrate(s)

operations (list):
    List of dictionaries, each describing one operation. Always with respect to the parent simulation if applicable 

inp (dict):
    property/value pairs describing the simulation input

output (dict):
    property/value pairs output by the calculation


For custom types ATOMS, NANOCLUSTER, ADSORBATE and SUBSTRATE see below.



Collection: workflows


_id (int): 
    unique identifier

username (str):
    user who executed the workflow
creation_time (str):
    time of creation of the workflow
parameters (dict):
    workflow-specific parameters
name (str):
    custom name of workflow
workflow_type (str):
    custom type of workflow



Collection: machine_learning


_id (int): 
    unique identifier

workflow_id (int):
    ID of workflow which the machine learning run was part of

method (str):
    name of the ML method: krr, nn, ...

method_params (dict):
    Parameters of the method

descriptor (str):
    name of the descriptor: soap, mbtr, cm, ...

descriptor_params (dict):
    Parameters of the descriptor used

training_set (int[]):
    list of simulation IDs used for training

validation_set (int[]):
    list of simulation IDs used in validation. If empty, cross-validation was used.
test_set (int[]):
    list of simulation IDs used in testing. If empty, only validation was used
prediction_set (int[]):
    list of simulation IDs used for prediction.
metrics_training (dict):
    dictionary of (“metric name”: value) on training set
    key: string = name of the metric
    value: float = calculated value
metrics_validation (dict):
    dictionary of (“metric name”: value) on validation set
    key: string = name of the metric
    value: float = calculated value
metrics_test (dict):
    dictionary of (“metric name”: value) on test set
    key: string = name of the metric
    value: float = calculated value
output (dict):
    relevant training output info


Ideally, method name corresponds to a python class/function in the platform, that is initialised with the parameter dictionary given in method_params. Similarly, descriptor name also matches a python class, to be initialised with its own given set of parameters, descriptor_params.
The field output is a dictionary with all the useful output values from the calculation.


Custom Type: ATOMS

A dictionary for describing atoms in a system, conceptually 
close to ase.Atoms object:


numbers (int[]):
    list of atomic numbers as numpy array [N] of ints
positions (float[N,3]):
    positions as numpy matrix [Nx3] of doubles
constraints (int[N,3]):
    frozen flags a matrix [Nx3] of int [optional] 1 = frozen, 0 = free
pbc (bool):
    use periodic boundaries
cell (float[3,3]):
    matrix 3x3 with cell vectors on the rows
celldisp (float[3,1]):
    displacement of cell from origin
info (dict):
    field for additional information related to structure


The order of atoms in this dictionary is the one found in the simulation input file.



Custom Type: ADSORBATE


reference_id (int):
    ID of the simulation to use as reference
atom_ids (int[]):
    atom indices in the ATOMS dictionary of the simulation record
site_class (str):
    class of adsorption site: “top”, “bridge”, “hollow”, “4-fold hollow”
site_ids (int[]):
    list of atom ids (in simulation record) that define the adsorption site


Custom Type: NANOCLUSTER

In general, simulation.nanoclusters is a list of dictionaries with this structure. 

reference_id (int):
    ID of the simulation where this cluster was made, -1 if original
atom_ids (int[]):
    atom indices in the ATOMS dictionary of the simulation record


Custom Type: SUBSTRATE


reference_id (int):
    ID of the parent support simulation, -1 if no parent
atom_ids (int[]):
    atom indices in the corresponding ATOMS dictionary. See below




Database query examples
-----------------------

A few examples how to query that database
are given in the gui/ folder on the github
repository.
