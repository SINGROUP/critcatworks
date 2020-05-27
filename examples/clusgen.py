from fireworks import LaunchPad, Workflow
import pathlib
import os,time, sys
import logging
import ase
import getpass

# internal modules
from critcatworks.workflows import generate_nanoclusters_workflow
from critcatworks.database import mylaunchpad

if __name__ == "__main__":
    IS_QUEUE = True
    USERNAME = "myusername"
    PASSWORD = getpass.getpass()
    if IS_QUEUE:
        logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logdir = str(pathlib.Path(".").resolve())
        logging.basicConfig(filename = logdir + "/nanocluster_workflow.log", level=logging.INFO)

    # set up the LaunchPad and reset it
    launchpad = mylaunchpad.create_launchpad(USERNAME, PASSWORD,)
    bondlength_dct = { 
        "Fe" : 2.3,
        "Co" : 2.3,
        "Ni" : 2.4,
        "Cu" : 2.5,
        "Ti" : 2.7,
        "Pt" : 2.6,}   
    wf = generate_nanoclusters_workflow(username = USERNAME, password = PASSWORD,
        worker_target_path = "/wrk/jagermar/DONOTREMOVE/workflow_runs/generator/testruns", 
        shape = "ico", nanocluster_size = 3, 
        compositions = [6, 13, 28, 42, 49], 
        elements = ["Fe","Co", "Ni","Cu","Ti", "Pt"], generate_pure_nanoclusters = True,
        n_configurations = 10, n_initial_configurations = 100, 
        bondlength_dct = bondlength_dct,
        extdb_connect = {"db_name": "ncdb"},
        )    

    # store workflow 
    launchpad.add_wf(wf)
