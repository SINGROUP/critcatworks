from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time, pathlib, sys
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import ase, ase.io
import logging, datetime
from critcatworks.database.extdb import update_workflows_collection
import numpy as np



@explicit_serialize
class InitialTask(FiretaskBase):
    """ 
    Task to initialize workflow database  """

    _fw_name = 'InitialTask'
    required_params = ['username', 'parameters', 'name', 'workflow_type']
    optional_params = []

    def run_task(self, fw_spec):

        username = self["username"]
        parameters = self["parameters"]
        name = self["name"]
        workflow_type = self["workflow_type"]

        creation_time = str(datetime.datetime.now(tz=None))

        #
        workflow = update_workflows_collection(username, creation_time, parameters = parameters,
            name = name, workflow_type = workflow_type)

        update_spec = fw_spec
        update_spec["temp"] = {}
        update_spec["simulations"] = {}
        update_spec["workflow"] = workflow
        update_spec["machine_learning"] = {}

        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)

def initialize_workflow_data(username, parameters, name = "UNNAMED", workflow_type = "UNNAMED"):
    firetask1  = InitialTask(username = username, parameters = parameters, name = name, workflow_type = workflow_type)
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'InitialTask'},
             name = 'InitialWork')
    return fw
