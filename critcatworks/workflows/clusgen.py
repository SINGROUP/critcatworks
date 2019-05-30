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

