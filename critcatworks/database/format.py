from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time, pathlib, sys
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import ase, ase.io


def atoms_dict_to_ase(atoms_dict):
    cell = atoms_dict['cell']
    celldisp = atoms_dict['celldisp']
    constraints = atoms_dict['constraints']
    pbc = atoms_dict['pbc']
    numbers = atoms_dict['numbers']
    positions = atoms_dict['positions']
    info = atoms_dict.get("info", {})

    atoms = ase.Atoms(numbers=numbers,
        positions=positions,
        cell = cell,
        celldisp = celldisp,
        constraint = constraints,
        pbc = pbc,
        info = info)
    return atoms


def ase_to_atoms_dict(atoms):
    positions = atoms.get_positions().tolist()
    cell = atoms.get_cell().tolist()
    pbc = atoms.get_pbc().tolist()
    numbers = atoms.get_atomic_numbers().tolist()
    try:
        constraints = atoms.constraints.tolist()
    except:
        constraints = atoms.constraints
    celldisp = atoms.get_celldisp().tolist()
    info = atoms.info

    atoms_dict = {
        "positions" : positions,
        "cell" : cell,
        "pbc" : pbc,
        "numbers" : numbers,
        "constraints" : constraints,
        "celldisp" : celldisp,
        "info" : info,    
        }
    return atoms_dict