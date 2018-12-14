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
from critcatworks.database import read_structures 
from critcatworks.dft.cp2k import setup_dry_folders, setup_dry_cp2k

def run_cp2k_dry(source_path, template_path, target_path = None):
    """
    Workflow to test if cp2k runs without problems
    """
    # FireWork: Read nanocluster structures and initialise a database
    # object containing set information
    abspath = pathlib.Path(source_path).resolve()

    if target_path == None:
        target_path = os.getcwd()
    else:
        target_path = pathlib.Path(target_path).resolve()

    fw_read_structures = read_structures(abspath)

    # Firework: setup folders for DFT calculations
    fw_setup_folders = setup_dry_folders(target_path = target_path)


    # add above Fireworks with links
    workflow_list = [fw_read_structures, 
        fw_setup_folders,
        ]

    links_dict = {
            fw_read_structures : [fw_setup_folders],
            }

    # FireWork: setup, run and extract DFT calculation
    # (involves checking for errors in DFT and rerunning)
    fw_setup_cp2k = setup_dry_cp2k(template_path = template_path)
    workflow_list.append(fw_setup_cp2k)
    links_dict[fw_setup_folders] = [fw_setup_cp2k]

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
    launchpad = LaunchPad(logdir=".", strm_lvl='INFO')
    launchpad.reset('', require_password=False)

    wf =  run_cp2k_dry(
        source_path = str(pathlib.Path("../../tests/dummy_db/nc_structures/").resolve()),
        template_path = str(pathlib.Path("../../tests/dummy_db/templates/nothing.template").resolve()), 
        target_path = str(pathlib.Path("../../tests/dummy_db/output/").resolve()),
        )

    # store workflow and launch it locally, single shot
    launchpad.add_wf(wf)

    # excecute workflow
    if IS_QUEUE:
        abspath = str(pathlib.Path(".").resolve())
        dft = CommonAdapter(
            q_type="SLURM",
            queue="test",
            nodes= 1,
            ntasks= 2,
            walltime= '00:01:00',
            constraint='hsw',
            account= None,
            job_name= 'dfttestrun',
            pre_rocket= "module load cp2k-env/4.1-hsw",
            post_rocket= "current fashion: post-modern rocket after running dft",
            logdir= abspath,
            #rocket_launch= "rlaunch  singleshot --offline")
            rocket_launch= "rlaunch  singleshot")

        for i in range(0, 1000):
            launch_rocket_to_queue(launchpad, FWorker(), dft, 
                launcher_dir=abspath + "/fw_logs", create_launcher_dir=True, reserve=True)
            time.sleep(3)
    else:
        #launch_rocket(launchpad, FWorker())
        rapidfire(launchpad, FWorker())

