# cp2k setup
# cp2k run
# cp2k parse/analysis
from fireworks import Firework, FWorker, LaunchPad, Workflow
import os,time, re, glob, sys, subprocess
from fireworks import explicit_serialize, FiretaskBase, FWAction
import pathlib, logging
import pycp2k, cp2kparser
import ase, ase.io
from critcatworks.database import atoms_dict_to_ase

@explicit_serialize
class CP2KSetupTask(FiretaskBase):
    """ 
    Task to setup DFT calculations.

    Args:
        None
    """

    _fw_name = 'CP2KSetupTask'
    required_params = ['template', 'target_path', 'calc_id', "n_max_restarts"]
    optional_params = ['name']

    def run_task(self, fw_spec):
        
        logging.info("CP2KSetupTask not implemented yet")
        
        #prefix = self.get("name", "cp2k_run_id")      
        template = self["template"]
        target_path = self["target_path"]
        calc_id = self["calc_id"]
        #n_max_restarts = self["n_max_restarts"]
        template_path = "template.txt"
        with open(template_path, "w") as the_file:
            the_file.write(template)
        # read template
        cp2kinput = glob.glob(template_path)[0]
        #calc = pycp2k.cp2k.CP2K()
        calc = pycp2k.CP2K()
        #inpparser = pycp2k.CP2KInputParser()
        #calc = inpparser.parse(calc, cp2kinput)
        calc.parse(cp2kinput)

        #logging.info("info about input parser")
        #logging.info(inpparser)
        #logging.info("cp2k info storage \n")
        #logging.info(inpparser.storage_obj)
        logging.debug("target_path")
        logging.debug(target_path)


        # TODO
        # get cell size of nanocluster / system
        # atoms = fw_spec["temp"]["calc_structures"]["calc_id"] 
        


        calc.CP2K_INPUT.FORCE_EVAL_list[0].SUBSYS.CELL.Abc = "[angstrom] 25 25 25"
        calc.CP2K_INPUT.FORCE_EVAL_list[0].SUBSYS.TOPOLOGY.Coord_file_name = "structure.xyz"
        calc.working_directory = str(target_path)
        logging.debug("working_directory: " + str(calc.working_directory))
        calc.project_name = "gopt"
        calc.write_input_file()
        logging.info("cp2k input file written TO" + calc.project_name + ".inp")


        input_string = calc.get_input_string()
        update_spec = fw_spec
        update_spec["input_string"] = input_string
        #pass_spec = fw_spec
        #print("dummy outputs generated")
        #fw_spec.pop("_category")
        #fw_spec.pop("name")
        #detours = Firework([CP2KRunTask(target_path=target_path, calc_id = calc_id, n_max_restarts = n_max_restarts)], 
        #    spec = {'_category' : "dft", 'name' : 'CP2KRunTask', "n_restarts" : 0},
        #    name = 'CP2KRunWork')
        #return FWAction(update_spec = fw_spec, detours = detours)
        return FWAction(update_spec = fw_spec)

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
        logging.info("Running CP2K")
        # shell command construction
        input_file = glob.glob(target_path + "/" + "*inp")[0]
        input_file = os.path.basename(input_file)
        output_file = input_file.replace(".inp", ".out")
        
        # check for restart file
        restart_file_list = glob.glob(target_path + "/" + "*restart")
        if len(restart_file_list) == 1:
            restart_file = restart_file_list[0]
            input_file = restart_file 
            input_file = os.path.basename(input_file)
        elif len(restart_file_list) > 1:
            logging.warning("Found several .restart files. Taking first one.")
            restart_file = restart_file_list[0]
            input_file = restart_file 
            input_file = os.path.basename(input_file)
        else:
            # otherwise use inp
            pass
        
        cp2k_bin="cp2k.popt"
        run_command = "srun " + cp2k_bin  + " -o " + output_file + " -i " + input_file
        command_list = run_command.split()
        print("shell command:")
        print(command_list)
        # running CP2K with srun in shell
        with cd(target_path):
            print("going into directory", target_path)
            subprocess.call(command_list, shell = False)

        print("run done")

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
        print("calc_id", calc_id)
        # read output

        # preparse for correct termination, warnings and exceeded walltime
        # preparsing needed for parser to run without throwing error
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
                total_energy = frame_sequence_potential_energy[-1]
            else:
                output_state = "not_converged"
                logging.info("CP2K not converged")
                print("CP2K not converged")
                total_energy = None

        # restart
        if output_state == "no_output" or output_state == "incorrect_termination" or output_state == "not_converged":
            if fw_spec["n_restarts"] < n_max_restarts:
                fw_spec["n_restarts"] += 1
                detours = Firework([CP2KRunTask(target_path=target_path, calc_id = calc_id, n_max_restarts = n_max_restarts)], 
                    spec = {'_category' : "dft", 'name' : 'CP2KRunTask', 'n_restarts' : int(n_restarts) + 1 },
                    name = 'CP2KRunWork')
                return FWAction(update_spec = fw_spec, detours = detours)    

        # update data
        if output_state == "ok" or output_state == "not_converged":
            logging.info("total_energy: " + str(total_energy))
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
    
            relaxed_structure = ase.Atoms(symbols = atom_labels[-1], positions = atom_positions[-1])
    
            atoms_dict = relaxed_structure.__dict__
    
            result_dict = {
                "is_converged" : is_converged,
                "number_of_frames_in_sequence" : number_of_frames_in_sequence,
                "nwarnings" : nwarnings,
                "is_walltime_exceeded" : is_walltime_exceeded,
                "total_energy" : total_energy,
            }

        else:
            # in case of no_output or incorrect_termination
            # fill slots with empty placeholders
            atoms_dict = {}
            input_file = ""
            result_dict = {"is_converged" : False, "output_state" : output_state}

        ##
        # get source simulation
        source_simulation = fw_spec["simulation"]

        # update external database
        dct = source_simulation.copy()
        # calculation originated from this:
        dct["source_id"] = calc_id
        dct["atoms"] = atoms_dict ### !!!
        dct["operations"] = ["cp2k"]
        dct["inp"]["input_string"] = fw_spec["input_string"]
        dct["output"] = result_dict # might still be missing some output

        simulation = update_simulations_collection(dct)

        # update internal workflow data
        simulation_id = simulation["_id"]
        mod_spec =[
            {"_set" : {"simulations->" + str(simulation_id) : simulation }},
            ]

        # update temp workflow data
    
            # mod_spec =[
            #     {'_set' : {'adsorbate_energies_dict->' + str(calc_id) : float(total_energy)}},
            #     {'_set' : {'is_converged_dict->' + str(calc_id) : is_converged}},
            #     {'_set' : {'relaxed_structure_dict->' + str(calc_id): atoms_dict}},
            #     {'_set' : {'dft_result_dict->' + str(calc_id) : result_dict}},
            #     {'_set' : {'n_restarts' : fw_spec["n_restarts"]}},
            #     ]


        fw_spec.pop("_category")
        fw_spec.pop("name")
        return FWAction(update_spec = fw_spec, mod_spec=mod_spec)


def setup_cp2k(template, target_path, calc_id, name = "cp2k_run_id", n_max_restarts = 4):
    setup_task = CP2KSetupTask(template = template,
                    target_path = target_path,
                    calc_id = calc_id,
                    name = name,
                    n_max_restarts = n_max_restarts)

    run_task = CP2KRunTask(target_path = target_path, 
        calc_id = calc_id, 
        n_max_restarts = n_max_restarts)
    
    # currently not part of setup to ensure robustness when cp2k fails
    #analysis_task = CP2KAnalysisTask(target_path = target_path, 
    #    calc_id = calc_id, 
    #    n_max_restarts = n_max_restarts)

    
    fw = Firework([setup_task,run_task], spec = {'_category' : "dft", 'name' : 'CP2KWork', "n_restarts" : 0},
                     name = 'CP2KWork')
    return fw


def rerun_cp2k(target_path, calc_id, n_max_restarts, n_restarts):
    """
    Run and analyse
    """
    fw = Firework([CP2KRunTask(target_path=target_path, calc_id = calc_id, n_max_restarts = n_max_restarts)], 
        spec = {'_category' : "dft", 'name' : 'CP2KRerunWork', 'n_restarts' : int(n_restarts) + 1 },
        name = 'CP2KRerunWork')
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
