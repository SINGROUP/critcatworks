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
        calc_paths = []
        #iterating over available structures
        for idx, atoms_dict in enumerate(ads_structures):
            atoms = atoms_dict_to_ase(atoms_dict)
            structure_folder = prefix + '_' + str(idx)

            structure_folder_path = parent_folder_path + "/" + structure_folder
            calc_paths.append(structure_folder_path)

            if not os.path.exists(structure_folder_path):
                os.makedirs(structure_folder_path)

            ase.io.write(structure_folder_path + "/" + "adsorbate_structure.xyz", atoms)
        
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
    required_params = ['template_path', 'target_path', 'chunk_size', 'name', 'n_max_restarts']
    optional_params = []

    def run_task(self, fw_spec):
        template_path = self["template_path"]
        target_path = self["target_path"]
        chunk_size = self["chunk_size"]
        name = self["name"]
        n_max_restarts = self["n_max_restarts"]
        calc_paths = fw_spec["calc_paths"]
        fps_ranking = fw_spec["fps_ranking"]

        # define what chunk to run
        try:
            n_calcs_started = fw_spec["n_calcs_started"]
        except KeyError:
            logging.info("Starting Chunk")
            n_calcs_started = 0

        ranked_ids = fps_ranking[n_calcs_started : n_calcs_started+chunk_size]

        detours = []
        for ranked_id in ranked_ids:
            logging.info(ranked_id)
            target_path = calc_paths[int(ranked_id)]

            # create detour to setup cp2k calculation
            new_fw = Firework([CP2KSetupTask(template_path = template_path,
                target_path = target_path,
                ranked_id = ranked_id,
                name = name,
                n_max_restarts = n_max_restarts)], spec = {'_category' : "lightweight", 'name' : 'CP2KSetupTask'},
                name = 'CP2KSetupWork')
            detours.append(new_fw)


        update_spec = fw_spec
        update_spec["n_calcs_started"] = n_calcs_started + chunk_size
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
    required_params = ['template_path', 'target_path', 'ranked_id', "n_max_restarts"]
    optional_params = ['name']

    def run_task(self, fw_spec):
        
        logging.info("CP2KSetupTask not implemented yet")
        
        prefix = self.get("name", "cp2k_run_id")      
        template_path = self["template_path"]
        target_path = self["target_path"]
        ranked_id = self["ranked_id"]
        n_max_restarts = self["n_max_restarts"]

        # read template
        cp2kinput = glob.glob(template_path)[0]
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
        #pass_spec = fw_spec
        #print("dummy outputs generated")
        fw_spec.pop("_category")
        fw_spec.pop("name")
        detours = Firework([CP2KRunTask(target_path=target_path, ranked_id = ranked_id, n_max_restarts = n_max_restarts)], 
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
    required_params = ['target_path', 'ranked_id', "n_max_restarts"]
    optional_params = []

    def run_task(self, fw_spec):
        target_path = self["target_path"]
        ranked_id = self["ranked_id"]
        n_max_restarts = self["n_max_restarts"]
        n_restarts = fw_spec["n_restarts"]
        logging.info("Running CP2K not implemented yet. Creating dummy outputs")
        print("Running CP2K dry. dummy outputs")
        print("ranked id", ranked_id)
        # shell command construction
        input_file = glob.glob(target_path + "/" + "*inp")[0]
        output_file = input_file.replace(".inp", ".out")
       
        
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
        
        #output_file ="dryrun.stuff"
        #print("redirecting cp2k output for testing to", output_file)
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
        detours = Firework([CP2KAnalysisTask(target_path=target_path, ranked_id = ranked_id, n_max_restarts = n_max_restarts)], 
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
    required_params = ['target_path', 'ranked_id', 'n_max_restarts']
    optional_params = []

    def run_task(self, fw_spec):
        target_path = self["target_path"]
        ranked_id = self["ranked_id"]
        n_max_restarts = self["n_max_restarts"]
        n_restarts = fw_spec["n_restarts"]
        print("ranked id", ranked_id)
        # read output

        # preparse for correct termination, warnings and exceeded walltime
        output_state, preparse_results = additional_parse_stdout(str(target_path))

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
            print("cp2k output file type", type(cp2koutput))
            results = parser.parse(str(cp2koutput))
            print("cp2k parser ran successfully")
    
            atom_positions = results["atom_positions"]
            number_of_frames_in_sequence = results["number_of_frames_in_sequence"]
            is_converged = results["geometry_optimization_converged"]
            atom_labels = results["atom_labels"]
            frame_sequence_potential_energy = results["frame_sequence_potential_energy"]
    
            if is_converged:
                adsorbate_total_energy = frame_sequence_potential_energy[-1]
            else:
                output_state = "not_converged"
                logging.info("CP2K not converged")
                print("CP2K not converged")
                adsorbate_total_energy = None

        if output_state == "ok":
            logging.info("adsorbate_total_energy: " + str(adsorbate_total_energy))
            print("adsorbate_total_energy: " + str(adsorbate_total_energy))
    
            logging.debug("atom_labels")
            logging.debug(atom_labels.shape)
            logging.debug(atom_labels[-1])
            logging.debug("atom_positions")
            logging.debug(atom_positions.shape)
            logging.debug("frame_sequence_potential_energy")
            logging.debug(frame_sequence_potential_energy.shape)
            logging.debug("number_of_frames_in_sequence")
            logging.debug(number_of_frames_in_sequence)
            logging.debug("is_converged")
            logging.debug(is_converged)
    
            print("is_converged")
            print(is_converged)
            relaxed_structure = ase.Atoms(symbols = atom_labels[-1], positions = atom_positions[-1])
    
            atoms_dict = relaxed_structure.__dict__
    
            result_dict = {
                "is_converged" : is_converged,
                "number_of_frames_in_sequence" : number_of_frames_in_sequence,
                "nwarnings" : nwarnings,
                "is_walltime_exceeded" : is_walltime_exceeded,
                "adsorbate_total_energy" : adsorbate_total_energy,
            }
    
            
            fw_spec.pop("_category")
            fw_spec.pop("name")
    
            mod_spec =[
                {'_set' : {'adsorbate_energies_dict->' + str(ranked_id) : float(adsorbate_total_energy)}},
                {'_set' : {'is_converged_dict->' + str(ranked_id) : 1}},
                {'_set' : {'relaxed_structure_dict->' + str(ranked_id): atoms_dict}},
                {'_set' : {'dft_result_dict->' + str(ranked_id) : result_dict}},
                {'_set' : {'n_restarts' : fw_spec["n_restarts"]}},
                ]
            return FWAction(update_spec = fw_spec, mod_spec=mod_spec)
        else:
            detours = Firework([CP2KRunTask(target_path=target_path, ranked_id = ranked_id, n_max_restarts = n_max_restarts)], 
                spec = {'_category' : "dft", 'name' : 'CP2KRunTask', 'n_restarts' : int(n_restarts) + 1 },
                name = 'CP2KRunWork')
            #keep track of number of restarts
            fw_spec["n_restarts"] += 1
            if fw_spec["n_restarts"] <= n_max_restarts:
                # restart
                return FWAction(update_spec = fw_spec, detours = detours)     
            else:
                result_dict = {
                    "is_converged" : False,
                    "n_restarts" : fw_spec["n_restarts"],
                }
                fw_spec.pop("_category")
                fw_spec.pop("name")
        
                mod_spec =[
                    {'_set' : {'adsorbate_energies_dict->' + str(ranked_id) : None}},
                    {'_set' : {'is_converged_dict->' + str(ranked_id) : 0}},
                    {'_set' : {'relaxed_structure_dict->' + str(ranked_id): None}},
                    {'_set' : {'dft_result_dict->' + str(ranked_id) : result_dict}},
                    {'_set' : {'n_restarts' : fw_spec["n_restarts"]}},
                    ]
                # no more restarts
                return FWAction(update_spec = fw_spec, mod_spec=mod_spec)


def setup_cp2k(template_path, target_path, chunk_size, name = "cp2k_run_id", n_max_restarts = 4):
    firetask1  = ChunkCalculationsTask(
        template_path = template_path,
        target_path = target_path,
        chunk_size = chunk_size,
        name = name,
        n_max_restarts = n_max_restarts,
        )
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'ChunkCalculationsTask'},
                     name = 'ChunkCalculationsWork')
    return fw


def setup_folders(target_path, name = "cp2k_run_id",):
    firetask1  = StructureFolderTask(
        target_path = target_path,
        name = name)
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'StructureFolderTask'},
             name = 'StructureFolderWork')
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
                                   

