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
    Task to update database from converged chunk
    of calculations.
    """

    _fw_name = 'UpdateDataTask'
    required_params = ['chunk_size']
    optional_params = []

    def run_task(self, fw_spec):
        chunk_size = int(self["chunk_size"])
        n_calcs_started = int(fw_spec["n_calcs_started"])

        ranked_ids_chunk = fw_spec["fps_ranking"][n_calcs_started - chunk_size : n_calcs_started]

        print(ranked_ids_chunk)

        adsorbate_energies_list = fw_spec["adsorbate_energies_list"]
        adsorbate_energies_dict = fw_spec["adsorbate_energies_dict"]

        relaxed_structure_list = fw_spec["relaxed_structure_list"]
        relaxed_structure_dict = fw_spec["relaxed_structure_dict"]


        logging.debug(adsorbate_energies_list)
        logging.debug(relaxed_structure_list)

        logging.debug(adsorbate_energies_dict)
        logging.debug(relaxed_structure_dict)


        for ranked_id in ranked_ids_chunk:
            adsorbate_energies_list[int(ranked_id)] = adsorbate_energies_dict[str(ranked_id)]
            relaxed_structure_list[int(ranked_id)] = relaxed_structure_dict[str(ranked_id)]

        update_spec = fw_spec
        update_spec["relaxed_structure_list"] = adsorbate_energies_list
        update_spec["adsorbate_energies_list"] = relaxed_structure_list


        logging.debug(adsorbate_energies_list)
        logging.debug(relaxed_structure_list)

        return FWAction(update_spec=update_spec)


def update_converged_data(chunk_size):
    firetask1  = UpdateDataTask(chunk_size = chunk_size)
    fw = Firework([firetask1])
    return fw
