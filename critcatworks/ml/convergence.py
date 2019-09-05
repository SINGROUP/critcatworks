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
    Task to check the convergence of the database.
    If not converged, the workflow continues.

    Args:
        threshold (float) : If the convergence_criterion (default MAE of property) is below the given threshold, 
                            the workflow is defused early
        convergence_criterion (str) :   Type of machine learning criterion, based on which to stop
                                        the workflow. Defaults to mae (MAE)

    Returns:
        FWAction : Firework action, updates fw_spec, possible defuses the workflow
    """

    _fw_name = 'CheckConvergenceTask'
    required_params = ['threshold', 'convergence_criterion']
    optional_params = []

    def run_task(self, fw_spec):

        threshold = self["threshold"]
        convergence_criterion = self["convergence_criterion"]
        
        machine_learning_id = fw_spec["temp"].get("last_machine_learning_id", "NONE")
        machine_learning = fw_spec["machine_learning"].get(str(machine_learning_id), {})
        mae = machine_learning.get("metrics_test", {}).get(convergence_criterion, 1000 * threshold)

        if mae < threshold:
            logging.info("Database is converged")
            logging.info("exiting workflow")
            defuse_children = True

        else:
            defuse_children = False
            # start new chunk of calculations
            # already defined as children in the workflow
            logging.info("calculating next chunk. continue workflow as planned")

        update_spec = fw_spec
        update_spec.pop("_category")
        update_spec.pop("name")

        return FWAction(update_spec=update_spec, defuse_children=defuse_children)


def check_convergence(threshold, convergence_criterion = "mae"):
    """
    Checks the convergence of the database.
    If not converged, the workflow continues.

    Args:
        threshold (float) : If the convergence_criterion (default MAE of property) is below the given threshold, 
                            the workflow is defused early
        convergence_criterion (str) :   Type of machine learning criterion, based on which to stop
                                        the workflow. Defaults to mae (MAE)

    Returns:
        Firework : Firework CheckConvergenceWork
    """
    firetask1  = CheckConvergenceTask(threshold = threshold, convergence_criterion = convergence_criterion)
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'CheckConvergenceTask'},
             name = 'CheckConvWork')
    return fw

