from fireworks import LaunchPad, Workflow
import pathlib
import os,time

# internal modules
from critcatworks.structure import generate_nanoclusters
from critcatworks.database.format import ase_to_atoms_dict
from critcatworks.database.update import initialize_workflow_data

def generate_nanoclusters_workflow(username, password, 
    worker_target_path = None, extdb_connect = {},
    shape = "ico", nanocluster_size = 3, compositions = "", elements = [], generate_pure_nanoclusters = True,
    n_configurations = 10, n_initial_configurations = 100, bondlength_dct = {},
    ):
    """
    Workflow to generate binary nanoclusters of defined size and shape.
    For each binary element combination, for each composition, n_configurations
    maximally dissimilar structures are created and uploaded to the
    simulations collection of the mongodb database.

    For generating the structures, the cluster generator in the 
    python package cluskit is used.

    Args:
        username (str) :        user who executed the workflow
        password (str) :        password for user to upload to the database
        worker_target_path (str) :  absolute path on computing resource 
                                    directory needs to exist
        extdb_connect (dict):   dictionary containing the keys host,
                                username, password, authsource and db_name.
        shape (str) : determines shape of nanoclusters. 'ico', 'octa' and 'wulff' 
        nanocluster_size (int) : determines nanocluster size. Meaning depends on shape 
        compositions (list) : each element determines the amount of atoms of type 1. 
        elements (list) : elements (str) to iterate over
        generate_pure_nanoclusters (bool) : if set to True, also pure 
                                            nanoclusters are generated
        n_configurations (int) :    number of configurations per composition 
                                    (chosen based on maximally different 
                                    structural features)
        n_initial_configurations (int) : number of initial configurations per 
                                         composition to choose from (higher number
                                         will make the grid finer)
        bondlength_dct (dict) :     bond lengths to use for specific elements. 
                                    Some default bond lenghts are provided for
                                    common elements

    Returns:
        fireworks.Workflow : generate_nanoclusters Fireworks Workflow object
    """

    #FireWork: Initialize workflow with workflow_id from external database
    parameters = {
        "worker_target_path" : worker_target_path,
        "workflow_type" : "generate_nanoclusters",
        "shape" : shape,
        "nanocluster_size" : nanocluster_size,
        "compositions" : compositions,
        "elements" : elements,
        "generate_pure_nanoclusters" : generate_pure_nanoclusters,
        "n_configurations" : n_configurations,
        "n_initial_configurations" : n_initial_configurations,
        "bondlength_dct" : bondlength_dct,
        }
    fw_init = initialize_workflow_data(username, password, parameters, 
        name = "UNNAMED", workflow_type = "generate_nanoclusters",
        extdb_connect = extdb_connect)

    # Firework: generate nanoclusters
    fw_generate_nanoclusters = generate_nanoclusters(n_initial_configurations, n_configurations, 
            shape, nanocluster_size, 
            compositions, elements, 
            generate_pure_nanoclusters = generate_pure_nanoclusters, bondlength_dct = bondlength_dct)

    # add above Fireworks with links
    workflow_list = [fw_init,
        fw_generate_nanoclusters,
        ]

    links_dict = {fw_init : [fw_generate_nanoclusters]}

    wf = Workflow(workflow_list, links_dict)
    return wf

