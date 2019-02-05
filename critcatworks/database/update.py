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
class GatherPropertyTask(FiretaskBase):
    """ 
    Task to update database from chunk
    of calculations. Computes properties of systems.
    """

    _fw_name = 'GatherPropertyTask'
    required_params = ['chunk_size']
    optional_params = ['adsite_types']

    def run_task(self, fw_spec):
        chunk_size = int(self["chunk_size"])
        adsite_types = self["adsite_types"]
        n_calcs_started = int(fw_spec["n_calcs_started"])

        calc_ids_chunk = fw_spec["temp"]["calc_ids"][n_calcs_started - chunk_size : n_calcs_started]
        logging.info(calc_ids_chunk)

        simulations = fw_spec["simulations"]

        # compute reaction energy and store them as lists for ml

        n_calcs = len(calc_ids)
        reaction_energies_list = np.zeros(n_calcs)
        is_converged_list = np.zeros(n_calcs)

        for calc_id in calc_ids_chunk:
            simulation = simulations[str(calc_id)]

            structure = simulation["atoms"]
            # TODO add how adsorbate moved
            # get closest site classified

            is_converged_list[int(calc_id)] = simulation["output"]["is_converged"]

            #if is_converged_list[int(calc_id)] == True:
            
            # get current simulation total_energy
            simulation_total_energy = simulation["output"]["total_energy"]
            # iterate over
            # adsorbates
            adsorbates = simulation["adsorbates"]
            # nanoclusters
            nanoclusters = simulation["nanoclusters"]
            # substrates
            substrates = simulation["substrates"]

            component_types = [adsorbates, nanoclusters, substrates]
            reaction_energy = simulation_total_energy
            for components in component_types:
                for component in components:
                    reference_id = component["reference_id"]
                    try:
                        total_energy = simulations[str(reference_id)]["output"]["total_energy"]
                    except:
                        logging.warning("total_energy not found!")
                        total_energy = 0.0
                    reaction_energy -= total_energy
                    reaction_energies_list[int(calc_id)] = reaction_energy

        update_spec = fw_spec

        update_spec["temp"]["property"] = reaction_energies_list
        update_spec["temp"]["is_converged_list"] = is_converged_list 

        update_spec.pop("_category")
        update_spec.pop("name")
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



def update_converged_data(chunk_size, adsite_types = ["top", "bridge", "hollow"]):
    firetask1  = GatherPropertyTask(chunk_size = chunk_size, adsite_types = adsite_types)
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'GatherPropertyTask'},
             name = 'GatherPropertyWork')
    return fw

def initialize_workflow_data(username, parameters, name = "UNNAMED", workflow_type = "UNNAMED"):
    firetask1  = InitialTask(username = username, parameters = parameters, name = name, workflow_type = workflow_type)
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'InitialTask'},
             name = 'InitialWork')
    return fw