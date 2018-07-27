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
from critcatworks.clusgeo import get_adsites
from critcatworks.database import read_structures

def get_adsites_workflow(path, adsorbate_energy=0.0, adsorbate_name='H'):
    """
    Workflow to determine the adsorption sites and energies of a set of
    nanocluster structures using CP2K and Clusgeo
    """
    # FireWork: Read nanocluster structures and initialise a database
    # object containing set information
    abspath = pathlib.Path(path).resolve()
    fw_read_structures = read_structures(abspath)
    # FireWork: Determine adsites and add to database
    fw_get_adsites = get_adsites(adsorbate_energy=0.0, adsorbate_name='H')
    # FireWork: FPS ranking

    # FireWork: setup, run and extract DFT calculation
    # (involves checking for errors in DFT and rerunning)

    # FireWork: update database, 
    # (includes reading relaxed structure and energy)

    # FireWork: machine learning from database

    # FireWork: check if converged, give intermediary overview.
    # give summary when finished


    wf = Workflow([fw_read_structures, fw_get_adsites], links_dict = {fw_read_structures: [fw_get_adsites]})
    return wf




if __name__ == "__main__":
    # set up the LaunchPad and reset it
    launchpad = LaunchPad()
    launchpad.reset('', require_password=False)
    #launchpad = LaunchPad(host="myhost", port=12345, \
    #name="fireworks_testing_db", username="my_user", \
    #password="my_pass")

    #wf = get_adsites_workflow()
    #wf = dummy_workflow()
    #wf = cp2k_test_workflow()
    #wf = test_foreachtask_workflow()
    wf = get_adsites_workflow(path = "/l/programs/critcatworks/tests/dummy_db/nc_structures/")

    # store workflow and launch it locally, single shot
    launchpad.add_wf(wf)


    # excecute workflow
    IS_QUEUE = False
    if IS_QUEUE:
        launch_rocket_to_queue(launchpad, FWorker(), adapter, launcher_dir=mypath, reserve=True)
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



