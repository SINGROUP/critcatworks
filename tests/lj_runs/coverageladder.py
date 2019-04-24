from fireworks import LaunchPad, Workflow
import pathlib
import os,time, sys
import logging
import ase
from scipy.spatial.distance import pdist
import getpass

# internal modules
from critcatworks.database import mylaunchpad
from critcatworks.workflows import get_coverage_ladder_workflow

if __name__ == "__main__":
    import logging
    IS_QUEUE = True
    USERNAME = "mjcritcat"
    #PASSWORD = getpass.getpass()
    PASSWORD = "heterogeniuscatalysis"
    if IS_QUEUE:
        logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logdir = str(pathlib.Path(".").resolve())
        logging.basicConfig(filename = logdir + "/coverage_workflow.log", level=logging.INFO)

    # set up the LaunchPad and reset it
    launchpad = mylaunchpad.create_launchpad(USERNAME, PASSWORD, server = "atlas")
    launchpad.reset('', require_password=False)
    
    wf = get_coverage_ladder_workflow(username = "mjcritcat", 
        password = PASSWORD,
        template_path = str(pathlib.Path("../templates/cp2k_mm_energy.inp").resolve()),     
        #worker_target_path = "/scratch/work/jagerm1/workflow_runs/coverage/production/selected_ptni_clusters_topbridgehollow_nremaining80",
        worker_target_path = "/wrk/jagermar/DONOTREMOVE/workflow_runs/coverage/testruns/selected_ptni_structures3/",
        start_ids = [7510],
        reference_energy = -1.16195386047558 * 0.5,
        adsorbate_name = "H",
        max_iterations = 4,
        n_max_restarts = 1,
        skip_dft = False,
        bond_length = 1.5,
        initial_direction = 1,
        ranking_metric = "similarity",
        d = 4,
        l = 2, 
        k = 7, 
        extdb_connect = {"db_name": "ncdb"},
    )

    # store workflow on launchpad
    launchpad.add_wf(wf)
