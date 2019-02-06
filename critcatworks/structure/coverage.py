from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time, pathlib, sys
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import ase, ase.io
import cluskit

import numpy as np
import logging

from critcatworks.database import atoms_dict_to_ase, ase_to_atoms_dict

@explicit_serialize
class AdsorbateEliminationTask(FiretaskBase):
    """ 
    Task to determine ranking of adsorption site structures.

    """

    _fw_name = 'AdsorbateEliminationTask'
    required_params = ['adsorbate_name', 'bond_length']
    optional_params = []

    def run_task(self, fw_spec):
        logging.debug(fw_spec)
        adsorbate_name = self['adsorbate_name']
        coverage_structures = fw_spec["coverage_structures"]

        # divide adsorbate and nanocluster
        coverage_structures_eliminated = []
        for idx, atoms_dict in enumerate(coverage_structures):
            atoms = atoms_dict_to_ase(atoms_dict)

            symbol_list = atoms.get_chemical_symbols()
            is_adsorbate = np.isin(symbol_list, adsorbate_name)
            is_cluster = np.isin(symbol_list, adsorbate_name, invert = True)
            adsorbate_atoms = atoms[is_adsorbate]
            cluster_atoms = atoms[is_cluster]
            adsorbate_positions = adsorbate_atoms.get_positions()
            
            kept_positions = cluskit.utils.x2_to_x(adsorbate_positions, bondlength = self['bond_length'])
            n_kept = kept_positions.shape[0]
            kept_adsorbate_atoms = ase.Atoms(symbols=[adsorbate_name] * n_kept, positions=kept_positions)

            kept_atoms = cluster_atoms + kept_adsorbate_atoms
            kept_atoms_dict = ase_to_atoms_dict(kept_atoms)
            coverage_structures_eliminated.append(kept_atoms_dict)
            print(kept_atoms_dict)

        update_spec = fw_spec
        update_spec["coverage_structures"] = coverage_structures_eliminated
        update_spec.pop("_category")
        update_spec.pop("name")

        return FWAction(update_spec=update_spec)




def eliminate_pairs(adsorbate_name, bond_length):
    firetask1  = AdsorbateEliminationTask(adsorbate_name = adsorbate_name, bond_length = bond_length)
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'AdsorbateEliminationTask'},
             name = 'AdsorbateEliminationWork')
    return fw
