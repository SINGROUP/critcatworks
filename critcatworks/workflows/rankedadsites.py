from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import pathlib

# internal modules
from critcatworks.clusgeo import get_adsites, rank_adsites
from critcatworks.database import read_structures, update_converged_data
from critcatworks.dft import setup_folders, setup_cp2k
from critcatworks.ml import get_mae, check_convergence

def get_adsites_workflow(source_path, template, target_path = None, adsorbate_energy=0.0, 
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
        adsorbate_energy=0.0, 
        adsorbate_name='H', 
        #adsite_types = ["top", "bridge", "hollow"],
        adsite_types = ["top"],
        )
    # FireWork: FPS ranking
    fw_rank_adsites = rank_adsites()

    # FireWork: setup, run and extract DFT calculation
    # (involves checking for errors in DFT and rerunning)
    fw_setup_folders = setup_folders(target_path = target_path)


    ### loop starts ###
    fw_setup_cp2k = setup_cp2k(template = template, target_path = target_path, chunk_size = chunk_size)
    # FireWork: update database, 
    # (includes reading relaxed structure and energy)

    fw_update_converged_data = update_converged_data(chunk_size = chunk_size)

    # FireWork: machine learning from database

    fw_get_mae = get_mae()

    # FireWork: check if converged, give intermediary overview.
    # give summary when finished

    fw_check_convergence = check_convergence(threshold = 0.1)
    ### loop ends ###

    max_iterations = int(max_calculations / chunk_size)

    workflow_list = [fw_read_structures, 
        fw_get_adsites, 
        fw_rank_adsites, 
        fw_setup_folders,
        fw_setup_cp2k,
        ]

    links_dict = {
            fw_read_structures: [fw_get_adsites], 
            fw_get_adsites: [fw_rank_adsites],
            fw_rank_adsites : [fw_setup_folders],
            fw_setup_folders : [fw_setup_cp2k],
            }

    for i in range(max_iterations):
        if i == 0:
            fw_update_converged_data = update_converged_data(chunk_size = chunk_size)
            fw_get_mae = get_mae()
            fw_check_convergence = check_convergence(threshold = 0.1)

            workflow_list.extend([fw_update_converged_data, fw_get_mae, fw_check_convergence])
        else:
            fw_setup_cp2k = setup_cp2k(template = template, target_path = target_path, chunk_size = chunk_size)
            fw_update_converged_data = update_converged_data(chunk_size = chunk_size)
            fw_get_mae = get_mae()
            fw_check_convergence = check_convergence(threshold = 0.1)


            workflow_list.extend([fw_setup_cp2k ,fw_update_converged_data, fw_get_mae, fw_check_convergence])

        links_dict = {
            fw_read_structures: [fw_get_adsites], 
            fw_get_adsites: [fw_rank_adsites],
            fw_rank_adsites : [fw_setup_folders],
            fw_setup_folders : [fw_setup_cp2k],
            }



    """
    wf = Workflow([fw_read_structures, 
        fw_get_adsites, 
        fw_rank_adsites, 
        fw_setup_folders, 
        fw_setup_cp2k,
        fw_update_converged_data,
        fw_get_mae,
        fw_check_convergence,
        ], 
        links_dict = {
            fw_read_structures: [fw_get_adsites], 
            fw_get_adsites: [fw_rank_adsites],
            fw_rank_adsites : [fw_setup_folders],
            fw_setup_folders : [fw_setup_cp2k],
            fw_setup_cp2k : [fw_update_converged_data],
            fw_update_converged_data : [fw_get_mae],
            fw_get_mae : [fw_check_convergence],
            #fw_check_convergence : [fw_setup_cp2k]
            })
    """
    wf = Workflow(workflow_list, links_dict)
    return wf




if __name__ == "__main__":
    import logging
    #logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    logging.basicConfig(level=logging.INFO)

    # set up the LaunchPad and reset it
    launchpad = LaunchPad()
    launchpad.reset('', require_password=False)
    #launchpad = LaunchPad(host="myhost", port=12345, \
    #name="fireworks_testing_db", username="my_user", \
    #password="my_pass")

    wf = get_adsites_workflow(
        source_path = "/l/programs/critcatworks/tests/dummy_db/nc_structures/", 
        #template = '/l/programs/critcatworks/tests/dummy_db/templates/cu_mm_bulk.inp', 
        template = '/l/programs/critcatworks/tests/dummy_db/templates/cp2k_mm_energy.inp', 
        target_path = "/l/programs/critcatworks/tests/dummy_db/output/", 
        chunk_size = 10,
        max_calculations = 10,
        )

    # store workflow and launch it locally, single shot
    launchpad.add_wf(wf)


    # excecute workflow
    IS_QUEUE = False
    if IS_QUEUE:
        launch_rocket_to_queue(launchpad, FWorker(), adapter, launcher_dir='.', reserve=True)
    else:
        #launch_rocket(launchpad, FWorker())
        rapidfire(launchpad, FWorker())


    if IS_QUEUE:
        # recover offline fireworks
        for i in range(0,10):
            time.sleep(5)
            ids =launchpad.get_fw_ids()
            for idx in ids:
                launchpad.recover_offline(launch_id = idx)



