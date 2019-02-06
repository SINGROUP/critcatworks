from fireworks import LaunchPad, Workflow

import pathlib
import os,time

# internal modules
from critcatworks.database.update import initialize_workflow_data
from critcatworks.database.format import ase_to_atoms_dict
from critcatworks.database import start_from_structures, start_from_database, update_converged_data
from critcatworks.structure import get_adsites, rank_adsites
from critcatworks.dft import setup_folders, chunk_calculations

from critcatworks.database import update_converged_data
from critcatworks.ml import get_mae, check_convergence


def get_singlesites_workflow(template_path, worker_target_path = None, structures = None, extdb_ids = None,
        source_path  = None, reference_energy=0.0,
        adsorbate_name='H', chunk_size = 100, max_calculations = 10000, username = "unknown", 
        adsite_types = ["top", "bridge", "hollow"], threshold = 0.1, n_max_restarts = 1):
    """

    Workflow to determine the adsorption sites and energies of a set of
    nanocluster structures using CP2K and ClusKit
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
        "max_calculations" : max_calculations,
        "adsorbate_name" : adsorbate_name,
        "chunk_size" : chunk_size,
        "adsite_types" : adsite_types,
        "descriptor" : "soap",
        "descriptor_params" : {"nmax" : 9, "lmax" :6, "rcut" : 5.0, 
            "crossover" : True, "sparse" : False},
        "simulation_method" : "cp2k",
        "threshold" : threshold,
        "n_max_restarts" : n_max_restarts,
        "workflow_type" : "singlesites",
        }

    fw_init = initialize_workflow_data(username, parameters, name = "UNNAMED", workflow_type = "singlesites")

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


    # FireWork: Determine adsites and add to database
    fw_get_adsites = get_adsites(
        reference_energy= reference_energy, 
        adsorbate_name= adsorbate_name, 
        adsite_types = adsite_types,
        descriptor = "soap",
        descriptor_params = {"nmax" : 9, "lmax" :6, "rcut" : 5.0, 
            "crossover" : True, "sparse" : False},
        )
    # FireWork: FPS ranking
    fw_rank_adsites = rank_adsites()

    # Firework: setup folders for DFT calculations
    fw_setup_folders = setup_folders(target_path = worker_target_path, name = "cp2k_singlesites_id")


    # add above Fireworks with links
    workflow_list = [fw_init,
        fw_get_structures, 
        fw_get_adsites, 
        fw_rank_adsites, 
        fw_setup_folders,
        ]

    links_dict = {
            fw_init: [fw_get_structures],
            fw_get_structures: [fw_get_adsites], 
            fw_get_adsites: [fw_rank_adsites],
            fw_rank_adsites : [fw_setup_folders],
            }

    ### loop starts ###
    max_iterations = int(max_calculations / chunk_size)
    for i in range(max_iterations):

        # FireWork: setup, run and extract DFT calculation
        # (involves checking for errors in DFT and rerunning)

        fw_chunk_calculations = chunk_calculations(template = template, target_path = worker_target_path, 
            chunk_size = chunk_size, n_max_restarts = n_max_restarts, simulation_method = "cp2k")
        workflow_list.append(fw_chunk_calculations)
        if i == 0:
            links_dict[fw_setup_folders] = [fw_chunk_calculations]
        else: 
            links_dict[fw_check_convergence] = [fw_chunk_calculations]

        # FireWork: update database, 
        # (includes reading relaxed structure and energy)
        fw_update_converged_data = update_converged_data(chunk_size = chunk_size)
        workflow_list.append(fw_update_converged_data)
        links_dict[fw_chunk_calculations] =[fw_update_converged_data]


        # FireWork: machine learning from database
        fw_get_mae = get_mae(target_path = worker_target_path)
        workflow_list.append(fw_get_mae)
        links_dict[fw_update_converged_data] =[fw_get_mae]


        # FireWork: check if converged, give intermediary overview.
        # give summary when finished
        fw_check_convergence = check_convergence(threshold = threshold)
        workflow_list.append(fw_check_convergence)
        links_dict[fw_get_mae] =[fw_check_convergence]

    ### loop ends ###

    wf = Workflow(workflow_list, links_dict)
    return wf

