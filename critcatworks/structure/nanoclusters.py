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
from critcatworks.database.extdb import update_simulations_collection
from critcatworks.database.format import atoms_dict_to_ase, ase_to_atoms_dict
import json
import numpy as np

@explicit_serialize
class NCStabilityTask(FiretaskBase):
    """ 
    Task to compare the stability of clusters using cohesive or total energy.
    Cohesive energy requires atomic total energies.

    Args:

    """

    _fw_name = 'NCStabilityTask'
    required_params = []
    optional_params = ['atomic_energies']

    def run_task(self, fw_spec):
        atomic_energies = self.get("atomic_energies", {})
        workflow_id = fw_spec.get("workflow", {"_id" : -1 }).get("_id", -1)
        logging.debug(fw_spec)
        update_spec = fw_spec
        analysis_ids = fw_spec["temp"]["analysis_ids"]
        # analysis_ids becomes calc_ids
        calc_ids = analysis_ids
        simulations = fw_spec["simulations"]

        cohesive_energy_dct = {}
        cohesive_energy_lst = []
        for calc_id in calc_ids:
            simulation = simulations[str(calc_id)]
            # get total energy
            try:
                total_energy = simulation["output"]["total_energy"]
            except:
                logging.warning("simulation" + str(calc_id) + " did not seam to have total_energy in its output.")
                continue

            atomic_numbers = simulation["atoms"]["numbers"]
            cohesive_energy = total_energy

            # get chemical formula
            dct = ase.utils.formula._count_symbols(atomic_numbers)
            for symbol, occurences in dct.items():
                cohesive_energy = cohesive_energy - ( atomic_energies.get(str(symbol), 0.0) * occurences )

            formula_string = ase.utils.formula.formula_metal(atomic_numbers, empirical=False)
            if formula_string not in cohesive_energy_dct:
                cohesive_energy_dct[formula_string] = {}
            cohesive_energy_dct[formula_string][str(calc_id)] = cohesive_energy
            cohesive_energy_lst.append(cohesive_energy)


        sort_ids = np.argsort(cohesive_energy_lst)
        sorted_calc_ids = np.array(calc_ids)[sort_ids]

        logging.info("cohesive_energy_dct")
        logging.info("cohesive_energy_lst")
        logging.info(cohesive_energy_dct)
        logging.info(cohesive_energy_lst)

        # make file with information ?
        # cluster_stability ?
        with open('cohesive_energy.json', 'w') as outfile:
            json.dump(cohesive_energy_dct, outfile)

        # preferrably updating workflow data


        # fireworks
        update_spec = fw_spec
        update_spec["temp"]["sorted_calc_ids"] = sorted_calc_ids
        update_spec["temp"]["cohesive_energy_dct"] = cohesive_energy_dct
        update_spec["calc_ids"] = calc_ids

        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)

def compare_nanoclusters(atomic_energies = {}):
    firetask1  = NCStabilityTask(atomic_energies = atomic_energies)
    dct = {'_category' : "lightweight", 'name' : 'NCStabilityTask'}
    fw = Firework([firetask1], spec=dct,
             name = 'NCStabilityWork')
    return fw
