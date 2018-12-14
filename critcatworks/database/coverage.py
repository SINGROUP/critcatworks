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

@explicit_serialize
class UpdateDataTask(FiretaskBase):
    """ 
    Task to update database from nanocluster coverage
    calculations.
    """

    _fw_name = 'UpdateDataTask'
    required_params = []
    optional_params = []

    def run_task(self, fw_spec):
        coverage_energies_dict = fw_spec["coverage_energies_dict"]
        history_coverage_structures_dict = fw_spec["history_coverage_structures_dict"]
        history_coverage_energies_dict = fw_spec["history_coverage_energies_dict"]
        relaxed_coverage_dict = fw_spec["relaxed_coverage_dict"]
        
        print("relaxed_structure_dict")
        print(relaxed_coverage_dict)

        print(coverage_energies_dict)
        logging.debug(coverage_energies_dict)
        logging.debug(relaxed_coverage_dict)
        coverage_structures = fw_spec["coverage_structures"]

        relaxed_structures_list = []
        # iterating over available structures
        # updating history
        # updating coverage_structures, 
        # which will be used to read structure for next step
        for idx, unrelaxed_atoms_dict in enumerate(coverage_structures):
            atoms_dict = relaxed_coverage_dict[str(idx)]
            if atoms_dict:
                pass
            else:
                atoms_dict = unrelaxed_atoms_dict
            relaxed_structures_list.append(atoms_dict)
            history_coverage_structures_dict[str(idx)].append(atoms_dict)
            history_coverage_energies_dict[str(idx)].append(coverage_energies_dict[str(idx)])

        update_spec = fw_spec
        update_spec["coverage_structures"] = relaxed_structures_list
        update_spec["history_coverage_structures_dict"] = history_coverage_structures_dict
        update_spec["history_coverage_energies_dict"] = history_coverage_energies_dict
        update_spec.pop("_category")
        update_spec.pop("name")

        return FWAction(update_spec=update_spec)


def update_coverage_data():
    firetask1  = UpdateDataTask()
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'UpdateDataTask'},
             name = 'UpdateDataWork')
    return fw
