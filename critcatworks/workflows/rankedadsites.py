from fireworks import LaunchPad, Workflow

import pathlib
import os,time
import mylaunchpad

# internal modules
from critcatworks.clusgeo import get_adsites, rank_adsites
from critcatworks.database import read_structures, update_converged_data
from critcatworks.dft import setup_folders, setup_cp2k
from critcatworks.ml import get_mae, check_convergence

def get_adsites_workflow(source_path, template_path, target_path = None, reference_energy=0.0, 
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
        fw_setup_cp2k = setup_cp2k(template_path = template_path, target_path = target_path, chunk_size = chunk_size, n_max_restarts = 1)
        workflow_list.append(fw_setup_cp2k)
        if i == 0:
            links_dict[fw_setup_folders] = [fw_setup_cp2k]
        else: 
            links_dict[fw_check_convergence] = [fw_setup_cp2k]

        # FireWork: update database, 
        # (includes reading relaxed structure and energy)
        fw_update_converged_data = update_converged_data(chunk_size = chunk_size)
        workflow_list.append(fw_update_converged_data)
        links_dict[fw_setup_cp2k] =[fw_update_converged_data]


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


if __name__ == "__main__":
    import logging
    IS_QUEUE = True
    if IS_QUEUE:
        logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logdir = str(pathlib.Path(".").resolve())
        logging.basicConfig(filename = logdir + "/logfile_ranked_adsites.log", level=logging.INFO)

    # set up the LaunchPad and reset it
    launchpad = mylaunchpad.create_launchpad()
    launchpad.reset('', require_password=False)

    wf = get_adsites_workflow(
        source_path = str(pathlib.Path("../../tests/dummy_db/nc_structures/").resolve()),
        template_path = str(pathlib.Path("../../tests/dummy_db/templates/rankedadsites_cheap_gopt.inp").resolve()), 
        target_path = str(pathlib.Path("../../tests/dummy_db/output/").resolve()),
        reference_energy = -1.16195386047558 * 0.5,
        adsorbate_name = "H",
        chunk_size = 5,
        max_calculations = 15,
        )

    # store workflow 
    launchpad.add_wf(wf)
