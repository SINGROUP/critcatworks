from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time, re, glob, sys, subprocess
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import pathlib, logging
import pycp2k, cp2kparser
import ase, ase.io
import shutil

from critcatworks.database import atoms_dict_to_ase

@explicit_serialize
class CoverageStructureFolderTask(FiretaskBase):
    """ 
    Task to setup folders with nanoclusters covered with
    adsorbates.

    Args:
        None
    """

    _fw_name = 'CoverageStructureFolderTask'
    required_params = ['target_path', 'name']
    optional_params = []

    def run_task(self, fw_spec):
        target_path = self['target_path']
        prefix = self['name']
        parent_folder_name = 'cp2k_calculations'
        parent_folder_path = target_path + "/" + parent_folder_name

        if not os.path.exists(parent_folder_path):
            os.makedirs(parent_folder_path)

        cov_structures = fw_spec["coverage_structures"]
        calc_paths = []
        #iterating over available structures
        for idx, atoms_dict in enumerate(cov_structures):
            atoms = atoms_dict_to_ase(atoms_dict)
            structure_folder = prefix + '_nc_' + str(idx)

            structure_folder_path = parent_folder_path + "/" + structure_folder
            calc_paths.append(structure_folder_path)

            if not os.path.exists(structure_folder_path):
                os.makedirs(structure_folder_path)

            ase.io.write(structure_folder_path + "/" + "coverage_structure.xyz", atoms)
        
        update_spec = fw_spec


        update_spec["calc_paths"] = calc_paths
        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec = update_spec)


@explicit_serialize
class ChunkCalculationsTask(FiretaskBase):
    """ 
    Create Fireworks with new calculations to setup and run

    Args:
        None
    """

    _fw_name = 'ChunkCalculationsTask'
    required_params = ['template_path', 'target_path', 'name', 'n_max_restarts']
    optional_params = []

    def run_task(self, fw_spec):
        template_path = self["template_path"]
        target_path = self["target_path"]
        name = self["name"]
        n_max_restarts = self["n_max_restarts"]
        calc_paths = fw_spec["calc_paths"]

        detours = []
        for calc_id, target_path in enumerate(calc_paths):
            # create detour to setup cp2k calculation
            new_fw = Firework([CP2KSetupTask(template_path = template_path,
                target_path = target_path,
                calc_id = calc_id,
                name = name,
                n_max_restarts = n_max_restarts)], spec = {'_category' : "lightweight", 'name' : 'CP2KSetupTask'},
                name = 'CP2KSetupWork')
            detours.append(new_fw)


        update_spec = fw_spec
        update_spec.pop("_category")
        update_spec.pop("name")

        return FWAction(update_spec = update_spec, detours = detours)


@explicit_serialize
class CP2KSetupTask(FiretaskBase):
    """ 
    Task to setup DFT calculations.

    Args:
        None
    """

    _fw_name = 'CP2KSetupTask'
    required_params = ['template_path', 'target_path', 'calc_id', "n_max_restarts"]
    optional_params = ['name']

    def run_task(self, fw_spec):
        
        logging.info("CP2KSetupTask not implemented yet")
        
        prefix = self.get("name", "cp2k_run_id")      
        template_path = self["template_path"]
        target_path = self["target_path"]
        calc_id = self["calc_id"]
        n_max_restarts = self["n_max_restarts"]

        # read template
        cp2kinput = glob.glob(template_path + "/" + "*inp")[0]
        print("cp2kinput")
        print(cp2kinput)
        is_copy_template = True
        if is_copy_template == True:
            shutil.copy(cp2kinput, str(target_path))  
        else:
            calc = pycp2k.cp2k.CP2K()
            inpparser = pycp2k.CP2KInputParser()
            calc = inpparser.parse(calc, cp2kinput)
    
            logging.info("info about input parser")
            logging.info(inpparser)
            logging.info("cp2k info storage \n")
            logging.info(inpparser.storage_obj)
            logging.debug("target_path")
            logging.debug(target_path)
            
            calc.CP2K_INPUT.FORCE_EVAL_list[0].SUBSYS.CELL.Abc = "[angstrom] 20 20 20"
    
            calc.working_directory = str(target_path)
            logging.debug("working_directory: " + str(calc.working_directory))
            calc.project_name = "gopt"
            calc.write_input_file()
              
            logging.info("cp2k input file written TO" + calc.project_name + ".inp")

        fw_spec.pop("_category")
        fw_spec.pop("name")
        detours = Firework([CP2KRunTask(target_path=target_path, calc_id = calc_id, n_max_restarts = n_max_restarts)], 
            spec = {'_category' : "dft", 'name' : 'CP2KRunTask', "n_restarts" : 0},
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
    required_params = ['target_path', 'calc_id', "n_max_restarts"]
    optional_params = []

    def run_task(self, fw_spec):
        target_path = self["target_path"]
        calc_id = self["calc_id"]
        n_max_restarts = self["n_max_restarts"]
        n_restarts = fw_spec["n_restarts"]
        print("calc id", calc_id)
        # shell command construction
        input_file = glob.glob(target_path + "/" + "*inp")[0]
        output_file = input_file.replace(".inp", ".out")
        output_file = os.path.basename(output_file)
        
        # check for restart file
        restart_file_list = glob.glob(target_path + "/" + "*restart")
        if len(restart_file_list) == 1:
            restart_file = restart_file_list[0]
            input_file = restart_file 
        elif len(restart_file_list) > 1:
            logging.warning("Found several .restart files. Taking first one.")
            restart_file = restart_file_list[0]
            input_file = restart_file 
        else:
            # otherwise use inp
            pass
        
        input_file = os.path.basename(input_file)
        cp2k_bin="cp2k.popt"
        run_command = "srun " + cp2k_bin  + " -o " + output_file + " -i " + input_file
        command_list = run_command.split()
        print("shell command:")
        print(command_list)
        # running CP2K with srun in shell
        with cd(target_path):
            print("going into directory", target_path)
            subprocess.call(command_list, shell = False)


        fw_spec.pop("_category")
        fw_spec.pop("name")
        detours = Firework([CP2KAnalysisTask(target_path=target_path, calc_id = calc_id, n_max_restarts = n_max_restarts)], 
            spec = {'_category' : "lightweight", 'name' : 'CP2KAnalysisTask', "n_restarts" : n_restarts},
            name = 'CP2KAnalysisWork')
        return FWAction(update_spec = fw_spec, detours = detours)


@explicit_serialize
class CP2KAnalysisTask(FiretaskBase):
    """ 
    Task to analyse CP2K calculations.

    Args:
        None
    """

    _fw_name = 'CP2KAnalysisTask'
    required_params = ['target_path', 'calc_id', 'n_max_restarts']
    optional_params = []

    def run_task(self, fw_spec):
        target_path = self["target_path"]
        calc_id = self["calc_id"]
        n_max_restarts = self["n_max_restarts"]
        n_restarts = fw_spec["n_restarts"]
        print("calc id", calc_id)
        # read output

        # preparse for correct termination, warnings and exceeded walltime
        output_state, preparse_results = additional_parse_stdout(str(target_path))
        print("output_state", output_state)
        if output_state == "no_output":
            logging.warning("no output file available")
            print("no output file available")
        elif output_state == "incorrect_termination":
            logging.info("incorrect termination")
            print("incorrect termination")
        else:
            logging.info("parser confirmed that CP2K terminated correctly")
            print("parser confirmed that CP2K terminated correctly")
        ####
        if output_state == "no_output" or output_state == "incorrect_termination":
            pass
        else:
            nwarnings = preparse_results['nwarnings']
            is_walltime_exceeded = preparse_results['exceeded_walltime']

            # cp2k parser
            parser = cp2kparser.CP2KParser(default_units=["hartree"], log_level=logging.INFO)
            print("cp2k parser setup successfully")
            cp2koutput = glob.glob(target_path + "/" + "*out")[0]
            print("cp2k output file", cp2koutput)
            results = parser.parse(str(cp2koutput))
            print("cp2k parser ran successfully")
    
            atom_positions = results["atom_positions"]
            number_of_frames_in_sequence = results["number_of_frames_in_sequence"]
            try:
                is_converged = results["geometry_optimization_converged"]
            except KeyError:
                is_converged = False            
            atom_labels = results["atom_labels"]
            frame_sequence_potential_energy = results["frame_sequence_potential_energy"]
            
            if is_converged:
                coverage_total_energy = frame_sequence_potential_energy[-1]
            else:
                output_state = "not_converged"
                logging.info("CP2K not converged")
                print("CP2K not converged")
                coverage_total_energy = None

        if output_state == "ok" or output_state == "not_converged":
            logging.info("coverage_total_energy: " + str(coverage_total_energy))
            print("coverage_total_energy: " + str(coverage_total_energy))
    
            logging.debug("atom_labels")
            logging.debug(atom_labels.shape)
            print(atom_labels.shape)
            logging.debug(atom_labels[-1])
            logging.debug("atom_positions")
            logging.debug(atom_positions.shape)
            print(atom_positions.shape)
            logging.debug("frame_sequence_potential_energy")
            logging.debug(frame_sequence_potential_energy.shape)
            logging.debug("number_of_frames_in_sequence")
            logging.debug(number_of_frames_in_sequence)
            logging.debug("is_converged")
            logging.debug(is_converged)
            print("atomic positions")
            print(atom_positions) 
            print("first atomic positions")
            print(atom_positions[0]) 
            print("last atomic positions")
            print(atom_positions[-1]) 
            print("is_converged")
            print(is_converged)
            # multiplier to convert from parser internal to angstrom
            relaxed_structure = ase.Atoms(symbols = atom_labels[-1], positions = atom_positions[-1] * 10e9)
    
            atoms_dict = relaxed_structure.__dict__
            print("relaxed structure from parser")
            print(atoms_dict)
    
            result_dict = {
                "is_converged" : is_converged,
                "number_of_frames_in_sequence" : number_of_frames_in_sequence,
                "nwarnings" : nwarnings,
                "is_walltime_exceeded" : is_walltime_exceeded,
                "coverage_total_energy" : coverage_total_energy,
            }
    
            
            fw_spec.pop("_category")
            fw_spec.pop("name")
    
            mod_spec =[
                {'_set' : {'coverage_energies_dict->' + str(calc_id) : str(coverage_total_energy)}},
                {'_set' : {'relaxed_coverage_dict->' + str(calc_id): atoms_dict}},
                {'_set' : {'dft_result_dict->' + str(calc_id) : result_dict}},
                ]
            return FWAction(update_spec = fw_spec, mod_spec=mod_spec)
        else:
            detours = Firework([CP2KRunTask(target_path=target_path, calc_id = calc_id, n_max_restarts = n_max_restarts)], 
                spec = {'_category' : "dft", 'name' : 'CP2KRunTask'},
                name = 'CP2KRunWork')
            #keep track of number of restarts
            fw_spec["n_restarts"] += 1
            if fw_spec["n_restarts"] <= n_max_restarts:
                # restart
                print("calculation has to be restarted")
                return FWAction(update_spec = fw_spec, detours = detours)     
            else:
                result_dict = {
                    "is_converged" : False,
                    "n_restarts" : fw_spec["n_restarts"],
                }
                fw_spec.pop("_category")
                fw_spec.pop("name")
        
                mod_spec =[
                    {'_set' : {'coverage_energies_dict->' + str(calc_id) : None}},
                    {'_set' : {'relaxed_coverage_dict->' + str(calc_id): None}},
                    {'_set' : {'dft_result_dict->' + str(calc_id) : result_dict}},
                    ]
                # no more restarts
                print("no more restarts!")
                return FWAction(update_spec = fw_spec, mod_spec=mod_spec)


def setup_coverage_cp2k(template_path, target_path, name = "cp2k_run_id", n_max_restarts = 4):
    firetask1  = ChunkCalculationsTask(
        template_path = template_path,
        target_path = target_path,
        name = name,
        n_max_restarts = n_max_restarts,
        )
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'ChunkCalculationsTask'},
                     name = 'ChunkCalculationsWork')
    return fw


def setup_coverage_folders(target_path, name = "cp2k_coverage",):
    firetask1  = CoverageStructureFolderTask(
        target_path = target_path,
        name = name)
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'CoverageStructureFolderTask'},
             name = 'CoverageStructureFolderWork')
    return fw


def additional_parse_stdout(target_path):
    output_files = glob.glob(target_path + "/" + "*out")
    if len(output_files) == 0:
        logging.warning("Cp2k output file not retrieved")
        return "no_output", None
    elif len(output_files) > 1:
        logging.warning("WARNING, more than one output file available")
    cp2koutput = output_files[0]
    abs_fn = pathlib.Path(cp2koutput).resolve()
    abs_fn = str(abs_fn)
    result_dict = {'exceeded_walltime': False}
    with open(abs_fn, "r") as f:
        for line in f.readlines():
            #if line.startswith(' ENERGY| '):
            #    result_dict['energy'] = float(line.split()[8])
            if 'The number of warnings for this run is' in line:
                result_dict['nwarnings'] = int(line.split()[-1])
            if 'exceeded requested execution time' in line:
                result_dict['exceeded_walltime'] = True

    if 'nwarnings' not in result_dict:
        logging.warning("CP2K did not finish properly.")
        return "incorrect_termination", result_dict
    else:
        return "ok", result_dict

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
                                   

