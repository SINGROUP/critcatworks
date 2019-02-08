from fireworks import LaunchPad, Workflow
import pathlib
import os,time

# internal modules
from critcatworks.database import start_from_structures, start_from_database, update_converged_data, 
from critcatworks.structures import compare_nanoclusters
from critcatworks.dft import setup_folders, chunk_calculations
from critcatworks.database.format import ase_to_atoms_dict
from critcatworks.database.update import initialize_workflow_data

def get_nanoclusters_workflow(template_path, worker_target_path = None, structures = None, 
    extdb_ids = None, source_path = None, reference_energy=0.0, username = "unknown", skip_dft = False):
    """
    Workflow to relax the structure of a set of
    nanoclusters using CP2K
    The given paths have to be absolute paths on the location where Firework jobs will be run.
    """

    with open (template_path, "r") as f:
        template = f.read()

    #FireWork: Initialize workflow with workflow_id from external database
    parameters = {
        "template" : template,
        "template_path" : template_path,
        "worker_target_path" : worker_target_path,
        "extdb_ids" : extdb_ids,
        "source_path" : source_path,
        "reference_energy" : reference_energy,
        }
    fw_init = initialize_workflow_data(username, parameters, name = "UNNAMED", workflow_type = "relax_nanoclusters")

    # FireWork: Read nanocluster structures and initialise a database
    # object containing set information
    if structures != None:
        jsonified_structures = []
        for atoms in structures:
            atoms_dict = ase_to_atoms_dict(atoms)
            jsonified_structures.append(atoms_dict) 
        fw_get_structures = start_from_structures(jsonified_structures)
    elif extdb_ids != None:
        fw_get_structures = start_from_database(extdb_ids)
    elif source_path != None:
        fw_get_structures = read_structures(source_path)
    else:
        raise ValueError('structures, extdb_ids or source_path contain no entries!')

    # Firework: setup folders for DFT calculations
    fw_setup_folders = setup_folders(target_path = worker_target_path, name = "cp2k_nanoclusters_id")


    # FireWork: setup, run and extract DFT calculation
    # (involves checking for errors in DFT and rerunning)
    fw_chunk_calculations = chunk_calculations(template = template, target_path = worker_target_path, 
        n_max_restarts = 1, simulation_method = "cp2k", skip_dft = skip_dft)

    # FireWork: compare stability of nanoclusters. 
    #Computes cohesive energies if atomic_energies of all involved elements are given
    fw_compare_nanoclusters =  compare_nanoclusters(atomic_energies = {})

    # add above Fireworks with links
    workflow_list = [fw_init,
        fw_get_structures, 
        fw_setup_folders,
        fw_chunk_calculations,
        fw_compare_nanoclusters,
        ]

    links_dict = {
            fw_init: [fw_get_structures],
            fw_get_structures: [fw_setup_folders],
            fw_setup_folders : [fw_chunk_calculations],
            fw_chunk_calculations : [fw_compare_nanoclusters],
            }

    wf = Workflow(workflow_list, links_dict)

    return wf

