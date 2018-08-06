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
    calculator = atoms_dict['_calc']
    cell = atoms_dict['_cell']
    celldisp = atoms_dict['_celldisp']
    constraint = atoms_dict['_constraints']
    pbc = atoms_dict['_pbc']
    numbers = atoms_dict['arrays']['numbers']
    positions = atoms_dict['arrays']['positions']
    info = atoms_dict['info']

    atoms = ase.Atoms(numbers=numbers,
        positions=positions,
        calculator=calculator,
        cell = cell,
        celldisp = celldisp,
        constraint = constraint,
        pbc = pbc,
        info = info)

    return atoms