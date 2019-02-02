from fireworks import LaunchPad, Workflow
import pathlib
import os,time

# internal modules
from critcatworks.database import start_from_structures, start_from_database, update_converged_data
from critcatworks.dft import setup_folders, chunk_calculations

def get_nanoclusters_workflow(template_path, worker_target_path = None, structures = None, extdb_ids = None, source_path = None, reference_energy=0.0):
    """
    Workflow to relax the structure of a set of
    nanoclusters using CP2K
    The given paths have to be absolute paths on the location where Firework jobs will be run.
    """
    with open (template_path, "r") as f:
        template = f.read()

    # FireWork: Read nanocluster structures and initialise a database
    # object containing set information
    if structures != None:
        jsonified_structures = []
        for atoms in structures:
            atoms_dict = atoms.__dict__
            jsonified_structures.append(atoms_dict) 
        fw_get_structures = start_from_structures(jsonified_structures)
    elif extdb_ids != None:
        fw_get_structures = start_from_database(extdb_ids)
    elif source_path != None:
        fw_get_structures = read_structures(source_path)
    else:
        raise ValueError('structures, extdb_ids or source_path contain no entries!')

    # Firework: setup folders for DFT calculations
    fw_setup_folders = setup_folders(target_path = worker_target_path, name = "cp2k_run_id")


    # FireWork: setup, run and extract DFT calculation
    # (involves checking for errors in DFT and rerunning)
    fw_chunk_calculations = chunk_calculations(template = template, target_path = worker_target_path, 
        n_max_restarts = 1, simulation_method = "cp2k")

    # add above Fireworks with links
    workflow_list = [fw_get_structures, 
        fw_setup_folders,
        fw_chunk_calculations
        ]

    links_dict = {
            fw_get_structures: [fw_setup_folders],
            fw_setup_folders : [fw_chunk_calculations]
            }

    wf = Workflow(workflow_list, links_dict)
    return wf

