from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time, pathlib, sys
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import ase, ase.io
import logging
import numpy as np
from critcatworks.dft.cp2k import setup_cp2k

@explicit_serialize
class CheckConvergenceTask(FiretaskBase):
    """ 
    Task to check the convergence of the database
    If not converged, the workflow continues.
    """

    _fw_name = 'CheckConvergenceTask'
    required_params = []
    optional_params = []

    def run_task(self, fw_spec):

        threshold = self["threshold"]
        convergence_criterion = self["convergence_criterion"]
        
        machine_learning_id = fw_spec["temp"]["last_machine_learning_id"]
        machine_learning = fw_spec["machine_learning"][str(machine_learning_id)]
        mae = machine_learning["metrics_test"][convergence_criterion]

        if mae < threshold:
            logging.info("Database is converged")
            logging.info("exiting workflow")
            defuse_workflow = True

        else:
            defuse_workflow = False
            # start new chunk of calculations
            # already defined as children in the workflow
            logging.info("calculating next chunk. continue workflow as planned")

        update_spec = fw_spec
        update_spec.pop("_category")
        update_spec.pop("name")

        return FWAction(update_spec=update_spec, defuse_workflow=defuse_workflow)


def check_convergence(threshold, convergence_criterion = "mae"):
    firetask1  = CheckConvergenceTask(threshold = threshold, convergence_criterion = convergence_criterion)
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'CheckConvergenceTask'},
             name = 'CheckConvergenceWork')
    return fw

