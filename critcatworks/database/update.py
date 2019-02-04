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

@explicit_serialize
class UpdateDataTask(FiretaskBase):
    """ 
    Task to update database from converged chunk
    of calculations.
    """

    _fw_name = 'UpdateDataTask'
    required_params = ['chunk_size']
    optional_params = []

    def run_task(self, fw_spec):
        chunk_size = int(self["chunk_size"])
        n_calcs_started = int(fw_spec["n_calcs_started"])
        nc_energies = fw_spec["nc_energies"]
        reverse_connect_dict = fw_spec["reverse_connect_dict"]
        reference_energy = fw_spec["reference_energy"]


        ranked_ids_chunk = fw_spec["fps_ranking"][n_calcs_started - chunk_size : n_calcs_started]

        logging.info(ranked_ids_chunk)

        adsorbate_energies_list = fw_spec["adsorbate_energies_list"]
        adsorbate_energies_dict = fw_spec["adsorbate_energies_dict"]

        is_converged_list = fw_spec["is_converged_list"]
        is_converged_dict = fw_spec["is_converged_dict"]

        relaxed_structure_list = fw_spec["relaxed_structure_list"]
        relaxed_structure_dict = fw_spec["relaxed_structure_dict"]

        reaction_energies_list = fw_spec["reaction_energies_list"]

        logging.debug(adsorbate_energies_list)
        logging.debug(relaxed_structure_list)

        logging.debug(adsorbate_energies_dict)
        logging.debug(relaxed_structure_dict)


        for ranked_id in ranked_ids_chunk:
            adsorbate_energies_list[int(ranked_id)] = adsorbate_energies_dict[str(ranked_id)]
            relaxed_structure_list[int(ranked_id)] = relaxed_structure_dict[str(ranked_id)]
            is_converged_list[int(ranked_id)] = is_converged_dict[str(ranked_id)]

        logging.debug("reverse_connect_dict")
        logging.debug(reverse_connect_dict)

        logging.info("reverse_connect_dict")
        pp(reverse_connect_dict)
        for ranked_id in ranked_ids_chunk:
            logging.info(ranked_id, type(ranked_id))
            nc_id = reverse_connect_dict[str(ranked_id)]

            if is_converged_list[int(ranked_id)] == True:
                reaction_energies_list[int(ranked_id)] = float(adsorbate_energies_list[int(ranked_id)]) - float(reference_energy) - float(nc_energies[int(nc_id)])
            else:
                reaction_energies_list[int(ranked_id)] = 0.0



        update_spec = fw_spec
        update_spec["relaxed_structure_list"] = relaxed_structure_list
        update_spec["adsorbate_energies_list"] = adsorbate_energies_list
        update_spec["reaction_energies_list"] = reaction_energies_list
        update_spec["is_converged_list"] = is_converged_list 
        update_spec.pop("_category")
        update_spec.pop("name")

        logging.debug(adsorbate_energies_list)
        logging.debug(relaxed_structure_list)

        return FWAction(update_spec=update_spec)


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



def update_converged_data(chunk_size):
    firetask1  = UpdateDataTask(chunk_size = chunk_size)
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'UpdateDataTask'},
             name = 'UpdateDataWork')
    return fw

def initialize_workflow_data(username, parameters, name = "UNNAMED", workflow_type = "UNNAMED"):
    firetask1  = InitialTask(username = username, parameters = parameters, name = name, workflow_type = workflow_type)
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'InitialTask'},
             name = 'InitialWork')
    return fw