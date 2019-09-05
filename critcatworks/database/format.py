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


def join_cluster_adsorbate(cluster, adsorbate):
    """
    Helper function to merge the structures
    cluster and adsorbate while retaining information
    about the ids

    Args:
        cluster (ase.Atoms) : nanocluster structure
        adsorbate (ase.Atoms) : single adsorbate

    Returns:
        tuple : ase.Atoms object of merged structure, ids of the
                nanocluster, ids of the adsorbate
    """
    joint_atoms = cluster + adsorbate
    cluster_ids = list(range(len(cluster)))
    adsorbate_ids = list(range(len(cluster_ids), len(joint_atoms)))

    return joint_atoms, cluster_ids, adsorbate_ids

def adsorbate_pos_to_atoms_lst(adspos, adsorbate_name):
    """
    Helper function to turn positions for adsorbates into
    ase atoms objects while the species is defined by
    adsorbate_name
    Attention! Works with only one adsorbate atom.
    In the future, cluskit might generalize to return a 
    list of adsorbates already in ase format.
    
    Args:
        adspos (2D ndarray) : positions of the adsorbate atoms
        adsorbate_name (str) : chemical symbol of the adsorbate atoms

    Returns:
        list : ase.Atoms objects of single atoms at each position
    """
    atoms_lst = []
    ads_structures_dict = []
    for adsorbate in adspos:
        atoms = ase.Atoms(symbols=adsorbate_name, positions=adsorbate.reshape((1,3)))
        atoms_lst.append(atoms)
    return atoms_lst    

def split_nanocluster_and_adsorbates(simulation):
        """Helper function to split a given structure in a simulation
        document into two ase atoms objects containing only
        the nanocluster or the adsorbate atoms.

        Args:
            simulation (dict) :

        Returns:
            tuple : cluster atoms (ase.Atoms),
                    adsorbate atoms (ase.Atoms),
                    adsorption site ids (list), 
                    adsorption site classes (list), 
                    reference ids (list), 
                    adsorbate ids (list)                    
        """
        adsorbate_ids = []
        reference_ids = []
        site_class_list = []
        site_ids_list = []
        atoms_dict = simulation["atoms"]
        atoms = atoms_dict_to_ase(atoms_dict)

        for adsorbate in simulation["adsorbates"]:
            adsorbate_ids.extend(adsorbate["atom_ids"])
            reference_ids.append(adsorbate["reference_id"])
            site_class_list.append(adsorbate.get("site_class", ""))
            site_ids_list.append(adsorbate.get("site_ids_list", []))
        adsorbate_atoms = atoms[np.array(adsorbate_ids, dtype=int)]

        cluster_ids = []
        for nanocluster in simulation["nanoclusters"]:
            cluster_ids.extend(nanocluster["atom_ids"])
        cluster_atoms = atoms[np.array(cluster_ids, dtype=int)]

        return cluster_atoms, adsorbate_atoms, site_ids_list, site_class_list, reference_ids, adsorbate_ids


