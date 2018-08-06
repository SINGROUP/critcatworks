from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import pathlib
import pycp2k
import ase, ase.io

from critcatworks.database import atoms_dict_to_ase

@explicit_serialize
class StructureFolderTask(FiretaskBase):
    """ 
    Task to setup folders with adsorbate structures.

    Args:
        None
    """

    _fw_name = 'StructureFolderTask'
    required_params = ['target_path', 'name']
    optional_params = []

    def run_task(self, fw_spec):
        target_path = self['target_path']
        prefix = self['name']
        parent_folder_name = 'cp2k_calculations'
        parent_folder_path = target_path + "/" + parent_folder_name

        if not os.path.exists(parent_folder_path):
            os.makedirs(parent_folder_path)

        ads_structures = fw_spec["ads_structures"]

        #iterating over available structures
        for idx, atoms_dict in enumerate(ads_structures):
            atoms = atoms_dict_to_ase(atoms_dict)
            structure_folder = prefix + '_' + str(idx)

            structure_folder_path = parent_folder_path + "/" + structure_folder
            

            if not os.path.exists(structure_folder_path):
                os.makedirs(structure_folder_path)

            ase.io.write(structure_folder_path + "/" + "adsorbate_structure.xyz", atoms)
            break

        return 

@explicit_serialize
class ChunkCalculationsTask(FiretaskBase):
    """ 
    Create Fireworks with new calculations to setup and run

    Args:
        None
    """

    _fw_name = 'ChunkCalculationsTask'
    required_params = ['template', 'target_path', 'chunk_size', 'name']
    optional_params = []

    def run_task(self, fw_spec):
        template = self["template"]
        target_path = self["target_path"]
        chunk_size = self["chunk_size"]
        name = self["name"]


        new_fws = CP2KSetupTask(template = template,
        target_path = target_path,
        chunk_size = chunk_size,
        name = name)

        return 


@explicit_serialize
class CP2KSetupTask(FiretaskBase):
    """ 
    Task to setup DFT calculations.

    Args:
        None
    """

    _fw_name = 'CP2KSetupTask'
    required_params = ['template', 'target_path']
    optional_params = ['name']

    def run_task(self, fw_spec):
        print("CP2KSetupTask not implemented yet")
        prefix = self.get("name", "cp2k_run_id")      
        template_path = self["template"]
        target_path = self["target_path"]

        # read template

        inpparser = pycp2k.CP2KInputParser()
        calc = inpparser.parse(template_path)

        print("info about input parser")
        print(inpparser)

        print("cp2k info storage \n")
        print(inpparser.storage_obj)



        #calc.working_directory = target_path
        #calc.project_name = "rewritten_cp2k"
        #calc.write_input_file()

        #print("file written TO", calc.project_name + ".inp")
        #pass_spec = fw_spec
        #print("dummy outputs generated")

        #print("TODO getting structures and energies")

        #mod_spec =[{'_set' : {'dft_outputs->id' + str(dft_params[0]) : dft_params[0]}}]
        return #FWAction(update_spec =pass_spec, mod_spec=mod_spec)



@explicit_serialize
class CP2KRunTask(FiretaskBase):
    """ 
    Task to run CP2K calculations.

    Args:
        None
    """

    _fw_name = 'CP2KRunTask'
    required_params = ['target_path']
    optional_params = []

    def run_task(self, fw_spec):
        print("Running CP2K not implemented yet")
        return

@explicit_serialize
class CP2KAnalysisTask(FiretaskBase):
    """ 
    Task to analyse CP2K calculations.

    Args:
        None
    """

    _fw_name = 'CP2KAnalysisTask'
    required_params = ['target_path']
    optional_params = []

    def run_task(self, fw_spec):
        print("Analysis of CP2K output not implemented yet")
        return




def setup_cp2k(template, target_path, chunk_size, name = "cp2k_run_id",):
    firetask1  = ChunkCalculationsTask(
        template = template,
        target_path = target_path,
        chunk_size = chunk_size,
        name = name,
        )
    fw = Firework([firetask1])
    return fw


def setup_folders(target_path, name = "cp2k_run_id",):
    firetask1  = StructureFolderTask(
        target_path = target_path,
        name = name)
    fw = Firework([firetask1])
    return fw




###############################################################################################

@explicit_serialize
class DFTSetupTask(FiretaskBase):
   """ 
   Task to setup DFT calculations.

   Args:
        None
   """

   _fw_name = 'DFTSetupTask'
   required_params = ['path']
   optional_params = []

   def run_task(self, fw_spec):
       #print(self.__dict__)
       print("DFTSetupTask not implemented yet")
       pp(fw_spec)
       print("forwarding information")
       pass_spec = fw_spec
       print("dummy outputs generated")

       print("TODO getting structures and energies")

       mod_spec =[{'_set' : {'dft_outputs->id' + str(dft_params[0]) : dft_params[0]}}]
       return FWAction(update_spec =pass_spec, mod_spec=mod_spec)



def branch_dft_calculations_fw():
    task = DFTSetupTask(spec={'dft_params':[]}, inputs=['dft_params'])
    task_dict = task.to_dict()
    firetask1 = ForeachTask(task=task_dict, split='dft_params', spec={'dft_params':[1,2,3]})
    fw = Firework([firetask1])
    return fw



def cp2k_test_workflow():
    simple_task = SimpleTestTask(spec={'simple': 'stuff'})
    firetask = PyTask(run_cp2k("cp2k_mm_energy.inp"))
    wf = Firework([simple_task, firetask])

    return wf


def run_cp2k(template_file):
    #export ASE_CP2K_COMMAND="mpirun -n 2 cp2k_shell.popt"
    from ase.calculators.cp2k import CP2K
    from ase.build import molecule

    with open('cp2k_mm_energy.inp', 'r') as f:
        content = f.read()
    calc = CP2K(inp=content)
    atoms = molecule('H2O', calculator=calc)
    atoms.center(vacuum=2.0)
    print(atoms.get_potential_energy())
    return
###############################################################################################