from fireworks import LaunchPad, Workflow

import pathlib
import os,time

# internal modules
from critcatworks.structure import start_coverage_ladder, add_remove_adsorbate, gather_ladder, step_coverage_ladder
from critcatworks.database.update import initialize_workflow_data
from critcatworks.database.format import ase_to_atoms_dict
from critcatworks.database import start_from_structures, start_from_database, read_structures
from critcatworks.structure import update_converged_data
from critcatworks.dft import setup_folders, chunk_calculations
from critcatworks.structure import get_per_type_coverage, eliminate_pairs, eliminate_closest


def get_coverage_ladder_workflow(template_path, username, password,
        worker_target_path = None, start_ids = None, reference_energy=0.0,
        adsorbate_name='H', max_iterations = 100, n_max_restarts = 1,
        skip_dft = False, bond_length = 1.5,
        d = 4, l = 2, k = 7, initial_direction = 1, ranking_metric = "similarity",
        extdb_connect = {}):
    """
    Workflow to determine a stable coverage of a nanocluster with single adsorbate atoms. One adsorbate at 
    a time is added or removed until certain break criteria are met.
    d, l, k, initial_direction and ranking_metric are parameters specific to the coverage ladder workflow.
    """
    
    with open (template_path, "r") as f:
        template = f.read()

    #FireWork: Initialize workflow with workflow_id from external database
    parameters = {
        "template" : template,
        "template_path" : template_path,
        "worker_target_path" : worker_target_path,
        "start_ids" : start_ids,
        "reference_energy" : reference_energy,
        "max_iterations" : max_iterations,
        "adsorbate_name" : adsorbate_name,
        "descriptor" : "soap",
        "descriptor_params" : {"nmax" : 9, "lmax" :6, "rcut" : 5.0, 
            "crossover" : True, "sparse" : False},
        "simulation_method" : "cp2k",
        "n_max_restarts" : n_max_restarts,
        "workflow_type" : "coverageladder",
        "d" : d,
        "l" : l,
        "k" : k,
        "bond_length" : bond_length,
        "ranking_metric" : ranking_metric,
        }

    fw_init = initialize_workflow_data(username, password, parameters, 
        name = "UNNAMED", workflow_type = "coverageladder",
        extdb_connect = extdb_connect)

    # FireWork: Initialize coverage ladder
    fw_start_coverage_ladder = start_coverage_ladder(start_ids, initial_direction = initial_direction)

    # add above Fireworks with links
    workflow_list = [fw_init,
        fw_start_coverage_ladder,
        ]

    links_dict = {
            fw_init : [fw_start_coverage_ladder],
            }

    ### loop starts ###
    for i in range(max_iterations):
        # Firework: add or remove one adsorbate, several times
        fw_add_remove_adsorbate =  add_remove_adsorbate(bond_length = bond_length, k = k, ranking_metric = ranking_metric)
        workflow_list.append(fw_add_remove_adsorbate)
        if i == 0:
            links_dict[fw_start_coverage_ladder] = [fw_add_remove_adsorbate]
        else:
            links_dict[fw_step_coverage_ladder] = [fw_add_remove_adsorbate]

        # Firework: setup folders for DFT calculations,
        fw_setup_folders = setup_folders(target_path = worker_target_path, name = "cp2k_ladder_iter_" + str(i))
        workflow_list.append(fw_setup_folders)
        links_dict[fw_add_remove_adsorbate] = [fw_setup_folders]

        # FireWork: setup, run and extract DFT calculation
        # (involves checking for errors in DFT and rerunning)
        fw_chunk_calculations = chunk_calculations(template = template, 
            target_path = worker_target_path, 
            chunk_size = -1, n_max_restarts = n_max_restarts, 
            simulation_method = "cp2k",
            name = "cp2k_ladder_iter_" + str(i), 
            skip_dft = skip_dft)
        workflow_list.append(fw_chunk_calculations)

        links_dict[fw_setup_folders] = [fw_chunk_calculations] 

        # FireWork: update database, 
        fw_gather_ladder = gather_ladder()
        workflow_list.append(fw_gather_ladder)
        links_dict[fw_chunk_calculations] = [fw_gather_ladder]

        # Firework: Decision to go up or down the coverage ladder
        fw_step_coverage_ladder = step_coverage_ladder(l = l, d = d,)
        workflow_list.append(fw_step_coverage_ladder)
        links_dict[fw_gather_ladder] =[fw_step_coverage_ladder]
    ### loop ends ###

    wf = Workflow(workflow_list, links_dict)
    return wf