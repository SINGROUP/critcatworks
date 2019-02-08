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
        calc_ids = fw_spec["temp"]["calc_ids"]
        # analysis ids becomes part of calc_ids
        analysis_ids = fw_spec["temp"]["analysis_ids"]
        simulations = fw_spec["simulations"]
        n_calcs = len(calc_ids)
        reaction_energies_list = fw_spec["temp"].get("property", np.zeros(n_calcs).tolist())
        is_converged_list = fw_spec["temp"].get("is_converged_list", np.zeros(n_calcs).tolist())

        calc_ids[n_calcs_started - chunk_size : n_calcs_started] = analysis_ids
        calc_ids_chunk = analysis_ids
        logging.info("Gather Properties of following calculations:")
        logging.info(calc_ids_chunk)

        # compute reaction energy and store them as lists for ml
        for idx, calc_id in zip(range(n_calcs_started - chunk_size, n_calcs_started), calc_ids_chunk):
            simulation = simulations[str(calc_id)]

            structure = simulation["atoms"]
            # TODO add how adsorbate moved
            # get closest site classified

            is_converged_list[idx] = simulation["output"]["is_converged"]

            print(is_converged_list[idx], idx)

            #if is_converged_list[int(calc_id)] == True:
            
            # get current simulation total_energy
            simulation_total_energy = simulation["output"].get("total_energy", 0.0)
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
                    reaction_energies_list[idx] = reaction_energy

        update_spec = fw_spec

        update_spec["temp"]["property"] = reaction_energies_list
        update_spec["temp"]["is_converged_list"] = is_converged_list 
        fw_spec["temp"]["analysis_ids"] = []
        fw_spec["temp"]["calc_ids"] = calc_ids
        print(is_converged_list)

        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)

def update_converged_data(chunk_size, adsite_types = ["top", "bridge", "hollow"]):
    firetask1  = GatherPropertyTask(chunk_size = chunk_size, adsite_types = adsite_types)
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'GatherPropertyTask'},
             name = 'GatherPropertyWork')
    return fw