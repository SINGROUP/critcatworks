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
class NCReadTask(FiretaskBase):
    """ 
    Task to setup DFT calculations.

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




def read_structures(path):
    firetask1  = NCReadTask(path=path)
    fw = Firework([firetask1], spec={'path':path, '_category' : "lightweight", 'name' : 'NCReadTask'},
             name = 'NCReadWork')
    return fw
