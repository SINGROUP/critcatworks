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
    Custom Firetask to initialize a new workflow instance 
    in the database.
    Additionally, initializes a few entries in the fw_spec.
    """

    _fw_name = 'InitialTask'
    required_params = ['username', 'password', 'parameters', 'name', 'workflow_type']
    optional_params = ['extdb_connect']

    def run_task(self, fw_spec):

        username = self["username"]
        password = self["password"]
        parameters = self["parameters"]
        extdb_connect = self["extdb_connect"]
        name = self["name"]
        workflow_type = self["workflow_type"]

        creation_time = str(datetime.datetime.now(tz=None))

        extdb_connect["username"] = username
        extdb_connect["password"] = password
        extdb_connect["host"] = extdb_connect.get("host",
            "nanolayers.dyndns.org:27017")

        extdb_connect["db_name"] = extdb_connect.get("db_name",
            "testdb")        
        extdb_connect["authsource"] = extdb_connect.get("authsource",
            extdb_connect["db_name"])

        workflow = update_workflows_collection(username, password, 
            creation_time, parameters = parameters,
            name = name, workflow_type = workflow_type, extdb_connect = extdb_connect)

        update_spec = fw_spec
        update_spec["temp"] = {}
        update_spec["simulations"] = {}
        update_spec["workflow"] = workflow
        update_spec["machine_learning"] = {}
        update_spec["extdb_connect"] = extdb_connect
        update_spec["temp"]["calc_analysis_ids_dict"] = {}
        update_spec["analysis_ids"] = []
        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)

def initialize_workflow_data(username, password, parameters, name = "UNNAMED", 
        workflow_type = "UNNAMED", extdb_connect = {}):
    """
    Creates a custom Firework object to initialize the workflow. 
    It updates the workflow collection and makes a few entries in 
    the fw_spec.

    Args:
        username (str) : username for the mongodb database
        password (str) : password for the mongodb database 
        parameters (dict) : workflow-specific input parameters
        name (str) :  custom name of the workflow
        workflow_type (str) :  custom workflow type
        extdb_connect (dict):   dictionary optionally containing the keys host,
                                authsource and db_name. All fields have a default
                                value.

    Returns:
        Firework object : InitialWork 
    """
    firetask1  = InitialTask(username = username, password = password, 
        parameters = parameters, name = name, 
        workflow_type = workflow_type, extdb_connect = extdb_connect)
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'InitialTask'},
             name = 'InitialWork')
    return fw
