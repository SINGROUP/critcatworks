from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

from pprint import pprint as pp
import pathlib
import os,time

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
        fw_setup_cp2k = setup_cp2k(template_path = template_path, target_path = target_path, chunk_size = chunk_size)
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
            fw_check_convergence : [fw_setup_cp2k]
            })
    """
    wf = Workflow(workflow_list, links_dict)
    return wf


if __name__ == "__main__":
    import logging
    #logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    logging.basicConfig(filename = "logfile_ranked_adsites.log", level=logging.INFO)

    # set up the LaunchPad and reset it
    launchpad = LaunchPad(logdir=".")
    launchpad.reset('', require_password=False)
    #launchpad = LaunchPad(host="myhost", port=12345, \
    #name="fireworks_testing_db", username="my_user", \
    #password="my_pass")

    wf = get_adsites_workflow(
        source_path = "/l/programs/critcatworks/tests/dummy_db/nc_structures/", 
        #template = '/l/programs/critcatworks/tests/dummy_db/templates/cu_mm_bulk.inp', 
        template_path = '/l/programs/critcatworks/tests/dummy_db/templates/', 
        target_path = "/l/programs/critcatworks/tests/dummy_db/output/",
        reference_energy = -1.16195386047558 * 0.5,
        adsorbate_name = "H",
        chunk_size = 12,
        max_calculations = 30,
        )

    # store workflow and launch it locally, single shot
    launchpad.add_wf(wf)


    # excecute workflow
    IS_QUEUE = False
    if IS_QUEUE:
        abspath = pathlib.Path(".").resolve()
        dft = CommonAdapter(
            q_type="SLURM",
            queue="test",
            nodes= 1,
            ntasks= 8,
            walltime= '00:02:00',
            constraint='hsw',
            account= None,
            job_name= 'dfttestrun',
            pre_rocket= None,
            post_rocket= None,
            logdir= abspath,
            rocket_launch= "rlaunch  singleshot --offline")
        lightweight = CommonAdapter(
            q_type="SLURM",
            queue="test",
            nodes= 1,
            ntasks= 8,
            walltime= '00:00:30',
            constraint='hsw',
            account= None,
            job_name= 'light',
            pre_rocket= None,
            post_rocket= None,
            logdir= abspath,
            rocket_launch= "rlaunch  singleshot --offline")
        medium = CommonAdapter(
            q_type="SLURM",
            queue="test",
            nodes= 1,
            ntasks= 8,
            walltime= '00:01:00',
            constraint='hsw',
            account= None,
            job_name= 'medium',
            pre_rocket= None,
            post_rocket= None,
            logdir= abspath,
            rocket_launch= "rlaunch  singleshot --offline")

        for i in range(0, 10):
            launch_rocket_to_queue(launchpad, FWorker(category='dft'), dft, launcher_dir=abspath, reserve=True)
            launch_rocket_to_queue(launchpad, FWorker(category='medium'), medium, launcher_dir=abspath, reserve=True)
            launch_rocket_to_queue(launchpad, FWorker(category='lightweight'), lightweight, launcher_dir=abspath, reserve=True)
    else:
        #launch_rocket(launchpad, FWorker())
        rapidfire(launchpad, FWorker(category=['dft', 'medium', 'lightweight']))
        #rapidfire(launchpad, FWorker(category='medium'))
        #rapidfire(launchpad, FWorker(category='lightweight'))


    # running in background to submit dynamic fireworks
    # and recover offline fireworks
    if IS_QUEUE:
        for i in range(0,10):
            # recover offline fireworks
            time.sleep(5)
            ids =launchpad.get_fw_ids()
            for idx in ids:
                launchpad.recover_offline(launch_id = idx)




