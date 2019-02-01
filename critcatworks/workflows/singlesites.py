from fireworks import LaunchPad, Workflow

import pathlib
import os,time
import mylaunchpad

# internal modules
from critcatworks.clusgeo import get_adsites, rank_adsites
from critcatworks.database import read_structures, update_converged_data
from critcatworks.dft import setup_folders, chunk_calculations
from critcatworks.ml import get_mae, check_convergence

def get_singlesites_workflow(source_path, template_path, target_path = None, reference_energy=0.0, 
        adsorbate_name='H', chunk_size = 100, max_calculations = 10000):
    """
    Workflow to determine the adsorption sites and energies of a set of
    nanocluster structures using CP2K and Clusgeo
    """
    # FireWork: Read nanocluster structures and initialise a database
    # object containing set information
    abspath = pathlib.Path(source_path).resolve()

    if target_path == None:
        target_path = os.getcwd()
    else:
        target_path = pathlib.Path(target_path).resolve()

    fw_read_structures = read_structures(abspath)
    # FireWork: Determine adsites and add to database
    fw_get_adsites = get_adsites(
        reference_energy=0.0, 
        adsorbate_name='H', 
        #adsite_types = ["top", "bridge", "hollow"],
        adsite_types = ["top"],
        )
    # FireWork: FPS ranking
    fw_rank_adsites = rank_adsites()

    # Firework: setup folders for DFT calculations
    fw_setup_folders = setup_folders(target_path = target_path)


    # add above Fireworks with links
    workflow_list = [fw_read_structures, 
        fw_get_adsites, 
        fw_rank_adsites, 
        fw_setup_folders,
        ]

    links_dict = {
            fw_read_structures: [fw_get_adsites], 
            fw_get_adsites: [fw_rank_adsites],
            fw_rank_adsites : [fw_setup_folders],
            }

    ### loop starts ###
    max_iterations = int(max_calculations / chunk_size)
    for i in range(max_iterations):

        # FireWork: setup, run and extract DFT calculation
        # (involves checking for errors in DFT and rerunning)
        fw_chunk_calculations = chunk_calculations(template_path = template_path, target_path = target_path, 
            chunk_size = chunk_size, n_max_restarts = 1, simulation_method = "cp2k")
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
        fw_get_mae = get_mae(target_path = target_path)
        workflow_list.append(fw_get_mae)
        links_dict[fw_update_converged_data] =[fw_get_mae]


        # FireWork: check if converged, give intermediary overview.
        # give summary when finished
        fw_check_convergence = check_convergence(threshold = 0.1)
        workflow_list.append(fw_check_convergence)
        links_dict[fw_get_mae] =[fw_check_convergence]

    ### loop ends ###

    wf = Workflow(workflow_list, links_dict)
    return wf

