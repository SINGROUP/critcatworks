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
from critcatworks.database.extdb import get_external_database
import numpy as np
from critcatworks.database.extdb import fetch_simulations

@explicit_serialize
class GatherPropertyTask(FiretaskBase):
    """ 
    Task to update database from chunk of calculations. 
    Computes properties of systems. Currently, only
    reaction energies (adsorption energies) are computed.
    The ids of the input structures in calc_ids are 
    replaced by the ids of the post-simulation structures 
    from analysis_ids.

    Args:
        chunk_size (int) :  number of calculations that are run simulataneously. 
                            Default -1 means all calculations are run at once.
        adsite_types (list) :   adsorption site types, can contain any combination of
                                "top", "bridge", "hollow"    

    Returns:
        FWAction : Firework action, updates fw_spec
    """
    _fw_name = 'GatherPropertyTask'
    required_params = ['chunk_size']
    optional_params = ['adsite_types']

    def run_task(self, fw_spec):
        calc_analysis_ids_dict = fw_spec["temp"]["calc_analysis_ids_dict"]
        chunk_size = int(self["chunk_size"])
        adsite_types = self["adsite_types"]
        n_calcs_started = int(fw_spec["n_calcs_started"])
        calc_ids = fw_spec["temp"]["calc_ids"]
        # analysis ids becomes part of calc_ids
        analysis_ids = fw_spec["temp"]["analysis_ids"]
        n_calcs = len(calc_ids)
        reaction_energies_list = fw_spec["temp"].get("property", np.zeros(n_calcs).tolist())
        is_converged_list = fw_spec["temp"].get("is_converged_list", np.zeros(n_calcs).tolist())

        # reorder analysis_ids
        reordered_analysis_ids = []
        for calc_id in calc_ids:
            if str(calc_id) in calc_analysis_ids_dict:
                analysis_id = calc_analysis_ids_dict[str(calc_id)]
                reordered_analysis_ids.append(analysis_id)

        analysis_ids =  reordered_analysis_ids

        print(chunk_size, type(chunk_size))
        if chunk_size == -1:
            calc_ids = analysis_ids
            id_range = range(len(calc_ids))
        else:
            calc_ids[n_calcs_started - chunk_size : n_calcs_started] = analysis_ids
            id_range = range(n_calcs_started - chunk_size, n_calcs_started)
        calc_ids_chunk = analysis_ids
        simulations = fetch_simulations(fw_spec["extdb_connect"], calc_ids_chunk)
        logging.info("Gather Properties of following calculations:")
        logging.info(calc_ids_chunk)

        ext_db =get_external_database(fw_spec["extdb_connect"])
        # compute reaction energy and store them as lists for ml
        print("id_range", id_range)
        for idx, calc_id in zip(id_range, calc_ids_chunk):
            simulation = simulations[str(calc_id)]

            structure = simulation["atoms"]
            # TODO add how adsorbate moved
            # get closest site classified

            is_converged_list[idx] = simulation["output"]["is_converged"]

            print(is_converged_list[idx], idx)
            
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
            print("energy before adding references", reaction_energy)
            for components in component_types:
                for component in components:
                    reference_id = component["reference_id"]
                    print(reference_id)
                    try:
                        reference_simulation = simulations[str(reference_id)]
                    except:
                        logging.info("getting reference from database")
                        reference_simulation  = ext_db["simulations"].find_one({"_id": reference_id})
                    try:
                        total_energy = reference_simulation["output"]["total_energy"]
                    except:
                        logging.warning("total_energy not found! Not contributing to reaction energy!")
                        total_energy = 0.0
                    try:
                        reaction_energy -= float(total_energy)
                    except:
                        logging.warning("Energy not understood!")
                        logging.warning(total_energy)

                    print(reaction_energy, "reference", reference_id)
            reaction_energies_list[idx] = reaction_energy

        fw_spec["temp"]["calc_analysis_ids_dict"] ={}
        fw_spec["temp"]["property"] = reaction_energies_list
        print("reaction_energies_list")
        print(reaction_energies_list)
        print(len(reaction_energies_list))
        fw_spec["temp"]["is_converged_list"] = is_converged_list 
        fw_spec["temp"]["analysis_ids"] = []
        fw_spec["temp"]["calc_ids"] = calc_ids
        print("is_converged_list")
        print(is_converged_list)
        print("calc_ids", calc_ids)

        fw_spec.pop("_category")
        fw_spec.pop("name")
        return FWAction(update_spec=fw_spec)

def update_converged_data(chunk_size, adsite_types = ["top", "bridge", "hollow"]):
    """ 
    Updates the internal fw_spec data from chunk of calculations. 
    Computes properties of systems. Currently, only
    reaction energies (adsorption energies) are computed.
    The ids of the input structures in calc_ids are 
    replaced by the ids of the post-simulation structures 
    from analysis_ids.

    Args:
        chunk_size (int) :  number of calculations that are run simulataneously. 
                            Default -1 means all calculations are run at once.
        adsite_types (list) :   adsorption site types, can contain any combination of
                                "top", "bridge", "hollow"    

    Returns:
        Firework : Firework GatherPropertyWork
    """
    firetask1  = GatherPropertyTask(chunk_size = chunk_size, adsite_types = adsite_types)
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'GatherPropertyTask'},
             name = 'GatherPropertyWork')
    return fw
