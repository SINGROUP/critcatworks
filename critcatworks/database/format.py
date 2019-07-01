from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time, pathlib, sys
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import ase, ase.io
import numpy as np

def atoms_dict_to_ase(atoms_dict):
    """
    Helper function to convert a ATOMS dictionary into an 
    ase.Atoms object

    Args:
        atoms_dict (dict) : dictionary with information about the atoms.
                            should be in the following format

                            numbers (1D ndarray)   : list of atomic numbers as numpy array [N] of ints
                            positions (2D ndarray) : positions as numpy matrix [Nx3] of doubles
                            constraints (2D ndarray) : frozen flags a matrix [Nx3] of int [optional] 1 = frozen, 0 = free
                            pbc (bool)             : use periodic boundaries
                            cell (2D ndarray)      : matrix 3x3 with cell vectors on the rows
                            celldisp (1D ndarray)  : displacement of cell from origin
                            info (dict)            : field for additional information related to structure

    Returns:
        ase.Atoms : Corresponding ase.Atoms object
    """
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
    """
    Helper function to convert an 
    ase.Atoms object into its corresponding python dictionary

    Args:
        atoms (ase.Atoms) : ase.Atoms object

    Returns:
        dict : Corresponding python dictionary
    """
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

def write_descmatrix(descmatrix):
    """
    Helper function to write a descriptor matrix required for 
    machine learning. It is stored as a file, since large arrays
    make fireworks slow.

    Args:
        descmatrix (2D np.ndarray) : descriptor matrix with 
                                     M features x N datapoints

    Returns:
        str : absolute path to file
    """
    time_str = time.strftime("%Y-%m-%d-%H-%M")
    path = "descmatrix_" + time_str +  ".npy"
    path = os.path.abspath(path)
    np.save(path, descmatrix)
    return path

def read_descmatrix(fw_spec):
    """
    Helper function to read a descriptor matrix required for 
    machine learning. It is stored as a file, since large arrays
    make fireworks slow.

    Args:
        fw_spec (dict) : Only the key 'descmatrix' is read. It expects a string 
                         with the absolute path to file    
                         
    Returns:
        2D np.ndarray  :    descriptor matrix with 
                            M features x N datapoints
    """
    path = fw_spec["temp"]["descmatrix"]
    descmatrix = np.load(path)
    return descmatrix


