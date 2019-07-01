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
from critcatworks.database.extdb import update_simulations_collection, get_external_database, _query_id_counter_and_increment
from critcatworks.database.format import atoms_dict_to_ase, ase_to_atoms_dict
import json
import numpy as np

@explicit_serialize
class NCReadTask(FiretaskBase):
    """ 
    Task to read nanocluster structures from xyz files.

    Args:
        Path (str)  :   Absolute path to a directory containing
                        structures readable by ASE
        cell_factor (float) :   enlarges cell size to x times the diameter
                                diameter of the structure

    Returns:
        FWAction : Firework action, update fw_spec
    """
    _fw_name = 'NCReadTask'
    required_params = ['path']
    optional_params = ['cell_factor']

    def run_task(self, fw_spec):
        path = self["path"]
        cell_factor = self.get("cell_factor", 2.5)

        logging.debug(fw_spec)

        structures = read_structures_locally(path, cell_factor = cell_factor)

        calc_ids = []
        for atoms_dict in ase_atoms_lst:
            atoms = atoms_dict_to_ase(atoms_dict)
            dct = atoms.info
            total_energy = dct.get("E", dct.get("energy", dct.get("total_energy", dct.get("TotalEnergy" , dct.get("totalenergy", "UNKNOWN")))))
            if total_energy == "UNKNOWN":
                logging.warning("total energy of cluster not specified in structure file")
            nanocluster_atom_ids = list(range(len(atoms)))

            db = get_external_database(fw_spec["extdb_connect"])
            simulations = db['simulations']
            # request id counter
            simulation_id = _query_id_counter_and_increment('simulations', db)

            nanocluster = {"atom_ids" : nanocluster_atom_ids, "reference_id" : simulation_id}

            dct = {"_id" : simulation_id, "atoms" : atoms_dict, 
                "source_id" : -1, "workflow_id" : workflow_id, 
                "nanoclusters" : [nanocluster], "adsorbates" : [], "substrates" : [], 
                "operations" : [], "inp" : {}, "output" : {"total_energy" : total_energy},
                }
            # update external database
            simulations.insert_one(dct)
            calc_ids.append(simulation_id)

        # update temp workflow data
        update_spec = fw_spec
        update_spec["temp"]["calc_ids"] = calc_ids
        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)


@explicit_serialize
class NCStartFromStructuresTask(FiretaskBase):
    """ 
    Task to setup starting structures from nanoclusters (ASE atoms objects).

    Args:
        ase_atoms_lst (str)  :   list of ASE atoms objects in dictionary format
    Returns:
        FWAction : Firework action, update fw_spec
    """

    _fw_name = 'NCStartFromStructuresTask'
    required_params = ['ase_atoms_lst']

    def run_task(self, fw_spec):
        ase_atoms_lst = self["ase_atoms_lst"]
        workflow_id = fw_spec.get("workflow", {"_id" : -1 }).get("_id", -1)
        logging.debug(fw_spec)
        update_spec = fw_spec

        calc_ids = []
        for atoms_dict in ase_atoms_lst:
            atoms = atoms_dict_to_ase(atoms_dict)
            dct = atoms.info
            total_energy = dct.get("E", dct.get("energy", dct.get("total_energy", dct.get("TotalEnergy" , dct.get("totalenergy", "UNKNOWN")))))
            if total_energy == "UNKNOWN":
                logging.warning("total energy of cluster not specified in structure file")
            nanocluster_atom_ids = list(range(len(atoms)))

            db = get_external_database(fw_spec["extdb_connect"])
            simulations = db['simulations']
            # request id counter
            simulation_id = _query_id_counter_and_increment('simulations', db)

            nanocluster = {"atom_ids" : nanocluster_atom_ids, "reference_id" : simulation_id}

            dct = {"_id" : simulation_id, "atoms" : atoms_dict, 
                "source_id" : -1, "workflow_id" : workflow_id, 
                "nanoclusters" : [nanocluster], "adsorbates" : [], "substrates" : [], 
                "operations" : [], "inp" : {}, "output" : {"total_energy" : total_energy},
                }
            # update external database
            simulations.insert_one(dct)
            calc_ids.append(simulation_id)

        # update temp workflow data
        update_spec["temp"]["calc_ids"] = calc_ids
        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)


@explicit_serialize
class NCStartFromDatabaseTask(FiretaskBase):
    """ 
    Task to setup starting structures from nanoclusters (ASE atoms objects).

    Args:
        db_ids_lst (str)  :     list of simulation ids in external database
        ext_db (pymongo)  :     external database pymongo object. Defaults to using
                                extdb_connect (dictionary containing the keys host, 
                                username, password, authsource and db_name).

    Returns:
        FWAction : Firework action, update fw_spec
    """

    _fw_name = 'NCStartFromDatabaseTask'
    required_params = ['db_ids_lst', 'ext_db']
    optional_params = []

    def run_task(self, fw_spec):
        db_ids_lst = self["db_ids_lst"]
        ext_db = self.get("ext_db", None)
        if ext_db == None:
            ext_db =get_external_database(fw_spec["extdb_connect"])
        workflow_id = fw_spec.get("workflow", {"_id" : -1 }).get("_id", -1)
        logging.debug(fw_spec)
        update_spec = fw_spec

        calc_ids = []
        for db_id in db_ids_lst:
            simulation = ext_db["simulations"].find_one({"_id": db_id})
            simulation["nanoclusters"][0]["reference_id"] = db_id
            simulation["source_id"] = db_id
            simulation["workflow_id"] = workflow_id

            dct = update_simulations_collection(extdb_connect = fw_spec["extdb_connect"],
                **simulation)
            # update internal workflow data
            simulation_id = dct["_id"]

            ## old update simulation internally.
            # update_spec["simulations"][str(simulation_id)] = dct
            calc_ids.append(simulation_id)

        # update temp workflow data
        update_spec["temp"]["calc_ids"] = calc_ids
        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)


def _empty_fields_spec(spec):
    """
    Helper function to return a dictionary with empty dictionaries:
    temp, simulations, workflow and machine_learning will be populated
    later.

    Args:
        spec (dict) : fw_spec dictionary

    Returns:
        dict :  the same given spec, unless it was empty. 
                In that case a few keys point to empty dictionaries
    """
    if spec == {}:
        return {"temp" : {}, "simulations" : {}, "workflow" : {}, 
            "machine_learning" : {}}
    else:
        return spec


def read_structures(path, spec = {}, cell_factor = 2.5):
    """
    Sets up Firework to read nanocluster structures from
    structure files (e.g xyz)
    In the second line of the input file, it looks for the keywords: 
    E, energy, total_energy, TotalEnergy, totalenergy
    It stores the first value found in the field output.total_energy

    The structures are stored in individual documents of the simulation 
    collection.

    Args:
        path (str) :    absolute path to the directory where the 
                        structure files (e.g. xyz format) can be found.
        spec (dict) :   optional additional entries for the fw_spec
        cell_factor (float) :   enlarges cell size to x times the diameter
                                diameter of the structure

    Returns:
        Firework : NCReadWork Firework
    """
    firetask1  = NCReadTask(path=path, cell_factor = cell_factor)
    dct = {'_category' : "lightweight", 'name' : 'NCReadTask'}
    dct.update(_empty_fields_spec(spec))
    fw = Firework([firetask1], spec=dct,
             name = 'NCReadWork')
    return fw

def start_from_structures(ase_atoms_lst, spec = {}):
    """
    Sets up Firework to read nanocluster structures from
    ASE atoms objects.
    The structures are copied to new individual documents of the simulation 
    collection. References to the current workflow, the parent nanocluster 
    and the source are updated.
    
    Args:
        ase_atoms_lst (str)  :   list of ASE atoms objects in dictionary format
        spec (dict) :   optional additional entries for the fw_spec

    Returns:
        FWAction : Firework action, update fw_spec
    """
    firetask1  = NCStartFromStructuresTask(ase_atoms_lst=ase_atoms_lst)
    dct = {'_category' : "lightweight", 'name' : 'NCStartFromStructuresTask'}
    dct.update(_empty_fields_spec(spec))
    fw = Firework([firetask1], spec=dct,
             name = 'NCStartFromStructuresWork')
    return fw    

def start_from_database(db_ids_lst, ext_db = None, spec = {}):
    """
    Sets up Firework to retrieve nanocluster structures from the 
    simulation collection of the mongodb database.
    In atoms.info it looks for the keywords: 
    E, energy, total_energy, TotalEnergy, totalenergy
    It stores the first value found in the field output.total_energy

    The structures are stored in individual documents of the simulation 
    collection.
    
    """
    firetask1  = NCStartFromDatabaseTask(db_ids_lst=db_ids_lst, ext_db = ext_db)
    dct = {'_category' : "lightweight", 'name' : 'NCStartFromDatabaseTask'}
    dct.update(_empty_fields_spec(spec))
    fw = Firework([firetask1], spec=dct,
             name = 'NCStartFromDatabaseWork')
    return fw


def read_structures_locally(path, cell_factor = 2.5):
    """
    Helper function to read structures locally. Can be used within
    a firework or outside.

    Args:
        path (str) :    absolute path to the directory where the 
                        structure files (e.g. xyz format) can be found.

        cell_factor (float) :   enlarges cell size to x times the diameter
                                diameter of the structure

    Returns:
        list :  list of ase.Atoms objects with a manipulated cellsize field.
    """
    structures = []
    path = pathlib.Path(path).resolve()
    for idx, p in enumerate(pathlib.Path(path).iterdir()):
        if p.is_file():
            logging.debug("nanocluster path " + str(p) + " stem " + str(p.stem))
            try:
                atoms = ase.io.read(str(p))
                # set cell to 2.5 (default) the diameter
                pos = atoms.get_positions()
                pdist(pos)
                diameter = pdist(pos).max()
                mpl = cell_factor
                
                atoms.set_cell([diameter * mpl, diameter * mpl, diameter * mpl])
                structures.append(atoms)
                logging.debug(atoms)
            except ValueError:
                logging.warning("WARNING: file type not understood" + str(p) )
                continue
            except:
                logging.error("Unexpected error:", sys.exc_info()[0])
    return structures