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
    USERNAME = "mjcritcat"
    #PASSWORD = getpass.getpass()
    PASSWORD = "heterogeniuscatalysis"
    if IS_QUEUE:
        logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logdir = str(pathlib.Path(".").resolve())
        logging.basicConfig(filename = logdir + "/nanocluster_workflow.log", level=logging.INFO)

    # set up the LaunchPad and reset it
    launchpad = mylaunchpad.create_launchpad(USERNAME, PASSWORD, server = "atlas")
    launchpad.reset('', require_password=False)

    wf = generate_nanoclusters_workflow(username = USERNAME, password = PASSWORD,
        worker_target_path = "/l/to_delete_temp/generate_nanoclusters", extdb_connect = {},
        shape = "ico", nanocluster_size = 3, compositions = [28,], elements = ["Ti", "W"], generate_pure_nanoclusters = True,
        n_configurations = 5, n_initial_configurations = 30, bondlength_dct = {},
        )    

    # store workflow 
    launchpad.add_wf(wf)
