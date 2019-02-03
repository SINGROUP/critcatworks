from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time, re, glob, sys
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import pathlib, logging
import pycp2k, cp2kparser
import ase, ase.io
import subprocess

from critcatworks.database import atoms_dict_to_ase, ase_to_atoms_dict

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

        nc_structures = fw_spec["nc_structures"]
        #iterating over available structures
        atoms_dict = nc_structures[0]
        atoms = atoms_dict_to_ase(atoms_dict)
        structure_folder = prefix
        structure_folder_path = parent_folder_path + "/" + structure_folder
        target_path = structure_folder_path

        if not os.path.exists(structure_folder_path):
            os.makedirs(structure_folder_path)

        ase.io.write(structure_folder_path + "/" + "adsorbate_structure.xyz", atoms)
        update_spec = fw_spec
        update_spec["target_path"] = target_path
        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec = update_spec)


@explicit_serialize
class CP2KSetupTask(FiretaskBase):
    """ 
    Task to setup DFT calculations.

    Args:
        None
    """

    _fw_name = 'CP2KSetupTask'
    required_params = ['template_path']
    optional_params = ['name']

    def run_task(self, fw_spec):
        
        logging.info("CP2KSetupTask not implemented yet")
        
        prefix = self.get("name", "cp2k_dry_run")      
        template_path = self["template_path"]
        target_path = fw_spec["target_path"]

        # read template
        cp2kinput = template_path
        calc = pycp2k.cp2k.CP2K()
        inpparser = pycp2k.CP2KInputParser()
        calc = inpparser.parse(calc, cp2kinput)

        logging.debug("target_path")
        logging.debug(target_path)
        
        calc.CP2K_INPUT.FORCE_EVAL_list[0].SUBSYS.CELL.Abc = "[angstrom] 20 20 20"

        calc.working_directory = str(target_path)
        logging.debug("working_directory: " + str(calc.working_directory))
        calc.project_name = "nonerun"
        calc.write_input_file()

        logging.info("cp2k input file written TO" + calc.project_name + ".inp")
        #pass_spec = fw_spec
        #print("dummy outputs generated")
        fw_spec.pop("_category")
        fw_spec.pop("name")
        detours = Firework([CP2KRunTask(target_path=target_path)], 
            spec = {'_category' : "dft", 'name' : 'CP2KRunTask'},
            name = 'CP2KRunWork')
        return FWAction(update_spec = fw_spec, detours = detours)


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
        target_path = self["target_path"]
        logging.info("Dry running CP2K")
        print("Dry running CP2K")

        # shell command construction
        input_file="nonerun.inp"
        output_file="nonerun.out"
        cp2k_bin="cp2k.popt"
        run_command = "srun " + cp2k_bin  + " -o " + output_file + " -i " + input_file
        command_list = run_command.split()
        print("shell command:")
        print(command_list)
        # running CP2K with srun in shell
        with cd(target_path):
            print("going into directory", target_path)
            subprocess.call(command_list, shell = False)

        print("dry run done")  

        fw_spec.pop("_category")
        fw_spec.pop("name")
        return FWAction(update_spec = fw_spec)



def setup_dry_cp2k(template_path, name = "cp2k_dry_run",):
    firetask1  = CP2KSetupTask(
        template_path = template_path,
        name = name,
        )
    fw = Firework([firetask1], spec = {'_category' : "dft", 'name' : 'CP2KSetupTask'},
             name = 'CP2KSetupWork')
    return fw


def setup_dry_folders(target_path, name = "cp2k_dry_run",):
    firetask1  = StructureFolderTask(
        target_path = target_path,
        name = name)
    fw = Firework([firetask1], spec = {'_category' : "dft", 'name' : 'StructureFolderTask'},
             name = 'StructureFolderWork')
    return fw


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
