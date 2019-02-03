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
from critcatworks.database.format import atoms_dict_to_ase

@explicit_serialize
class NCReadTask(FiretaskBase):
    """ 
    Task to read nanocluster structures from xyz files.

    Args:
        Path (str)  :   Absolute path to a directory containing
                        structures readable by ASE
    """

    _fw_name = 'NCReadTask'
    required_params = ['path']
    optional_params = []

    def run_task(self, fw_spec):
        path = self["path"]
        logging.debug(fw_spec)

        f_lst = []
        nc_names = []
        nc_structures_ase = []
        nc_structures_dict = []
        nc_energies = []
        nc_ids = []
        atomic_numbers = []
        for idx, p in enumerate(pathlib.Path(path).iterdir()):
            if p.is_file():
                logging.debug("nanocluster path " + str(p) + " stem " + str(p.stem))
                f_lst.append(p)
                try:
                    atoms = ase.io.read(str(p))

                    logging.debug(atoms)
                except ValueError:
                    logging.warning("WARNING: file type not understood" + str(p) )
                    continue
                except:
                    logging.error("Unexpected error:", sys.exc_info()[0])
                nc_names.append(p.stem)
                atomic_numbers.extend(atoms.get_atomic_numbers())
                atoms_dict = atoms.__dict__
                nc_structures_ase.append(atoms)
                nc_structures_dict.append(atoms_dict)
                nc_ids.append(idx)

                energy = None
                energy = atoms.info.get("E")
                nc_energies.append(energy)

        sorted_set_atomic_numbers = sorted(set(atomic_numbers)) 
        # create dictionary entries
        update_spec = fw_spec
        update_spec["nc_structures"] = nc_structures_dict
        update_spec["nc_names"] = nc_names
        update_spec["nc_ids"] = nc_ids
        update_spec["nc_energies"] = nc_energies
        update_spec["nc_atomic_numbers"] = sorted_set_atomic_numbers
        update_spec.pop("_category")
        update_spec.pop("name")

        return FWAction(update_spec=update_spec)


@explicit_serialize
class NCStartFromStructuresTask(FiretaskBase):
    """ 
    Task to setup starting structures from nanoclusters (ASE atoms objects).

    Args:
        ase_atoms_lst (str)  :   list of ASE atoms objects in dictionary format
    """

    _fw_name = 'NCStartFromStructuresTask'
    required_params = ['ase_atoms_lst']
    optional_params = ['ext_db']

    def run_task(self, fw_spec):
        ase_atoms_lst = self["ase_atoms_lst"]
        ext_db = self.get("ext_db", None)
        workflow_id = fw_spec.get("workflow", {"_id" : -1 }).get("_id", -1)
        logging.debug(fw_spec)
        update_spec = fw_spec

        calc_ids = []
        for atoms_dict in ase_atoms_lst:
            atoms = atoms_dict_to_ase(atoms_dict)
            nanocluster_atom_ids = list(range(len(atoms)))
            nanocluster = {"atom_ids" : nanocluster_atom_ids, "reference_id" : -1}

            # update external database
            # enter datapoint
            dct = update_simulations_collection(wf_sim_id = -1, 
                atoms = atoms_dict, 
                source_id = -1, workflow_id = workflow_id, 
                nanoclusters = [nanocluster], adsorbates = [], substrates = [], 
                operations = [], inp = {}, output = {},
                )

            # update internal workflow data
            simulation_id = dct["_id"]

            # update simulation internally.
            dct["nanoclusters"][0]["reference_id"] = simulation_id

            update_spec["simulations"][str(simulation_id)] = dct
            calc_ids.append(simulation_id)

        # update temp workflow data
        update_spec["temp"]["calc_ids"] = calc_ids

        # collective data easier for consecutive processing

        # fireworks
        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)


@explicit_serialize
class NCStartFromDatabaseTask(FiretaskBase):
    """ 
    Task to setup starting structures from nanoclusters (ASE atoms objects).

    Args:
        db_ids_lst (str)  :     list of simulation ids in external database
        ext_db (pymongo)  :     external database pymongo object
    """

    _fw_name = 'NCStartFromDatabaseTask'
    required_params = ['db_ids_lst', 'ext_db']
    optional_params = []

    def run_task(self, fw_spec):
        db_ids_lst = self["db_ids_lst"]
        ext_db = self.get("ext_db", None)
        workflow_id = fw_spec.get("workflow", {"_id" : -1 }).get("_id", -1)
        logging.debug(fw_spec)
        update_spec = fw_spec

        calc_ids = []
        for db_id in db_ids_lst:
            simulation = ext_db["simulations"].find_one({"_id": db_id})
            
            # update internal workflow data
            simulation_id = simulation["_id"]

            update_spec["simulations"][str(simulation_id)] = simulation
            calc_ids.append(simulation_id)

        # update temp workflow data
        update_spec["temp"]["calc_ids"] = calc_ids
        # collective data easier for consecutive processing

        # fireworks
        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)

def _empty_fields_spec(spec):
    if spec == {}:
        return {"temp" : {}, "simulations" : {}, "workflow" : {}, 
            "machine_learning" : {}}
    else:
        return spec


def read_structures(path, spec = {}):
    firetask1  = NCReadTask(path=path)
    dct = {'_category' : "lightweight", 'name' : 'NCReadTask'}
    dct.update(_empty_fields_spec(spec))
    fw = Firework([firetask1], spec=dct,
             name = 'NCReadWork')
    return fw

def start_from_structures(ase_atoms_lst, ext_db = None, spec = {}):
    firetask1  = NCStartFromStructuresTask(ase_atoms_lst=ase_atoms_lst, ext_db = ext_db)
    dct = {'_category' : "lightweight", 'name' : 'NCStartFromStructuresTask'}
    dct.update(_empty_fields_spec(spec))
    fw = Firework([firetask1], spec=dct,
             name = 'NCStartFromStructuresWork')
    return fw    

def start_from_database(db_ids_lst, ext_db = None, spec = {}):
    firetask1  = NCStartFromDatabaseTask(db_ids_lst=db_ids_lst, ext_db = ext_db)
    dct = {'_category' : "lightweight", 'name' : 'NCStartFromDatabaseTask'}
    dct.update(_empty_fields_spec(spec))
    fw = Firework([firetask1], spec=dct,
             name = 'NCStartFromDatabaseWork')
    return fw    