from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time, pathlib, sys
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import ase, ase.io


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
        pp(fw_spec)

        f_lst = []
        nc_names = []
        nc_structures_ase = []
        nc_structures_dict = []
        nc_energies = []
        nc_ids = []
        for idx, p in enumerate(pathlib.Path(path).iterdir()):
            if p.is_file():
                print(p, p.stem)
                f_lst.append(p)
                try:
                    atoms = ase.io.read(p)

                    print(atoms)
                except ValueError:
                    print("WARNING: file type not understood", p)
                    continue
                except:
                    print("Unexpected error:", sys.exc_info()[0])
                nc_names.append(p.stem)
                atoms_dict = atoms.__dict__
                nc_structures_ase.append(atoms)
                nc_structures_dict.append(atoms_dict)
                nc_ids.append(idx)
                try:
                    energy = atoms.info["energy"]
                except:
                    energy = None

                nc_energies.append(energy)

        # create dictionary entries
        update_spec = fw_spec
        update_spec["nc_structures"] = nc_structures_dict
        update_spec["nc_names"] = nc_names
        update_spec["nc_ids"] = nc_ids
        update_spec["nc_energies"] = nc_energies

        return FWAction(update_spec=update_spec)




def read_structures(path):
    firetask1  = NCReadTask(path=path)
    fw = Firework([firetask1], spec={'path':path})
    return fw
