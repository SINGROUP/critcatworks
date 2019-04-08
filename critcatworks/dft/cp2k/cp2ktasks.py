# cp2k setup
# cp2k run
# cp2k parse/analysis
from fireworks import Firework, FWorker, LaunchPad, Workflow
import os,time, re, glob, sys, subprocess
from fireworks import explicit_serialize, FiretaskBase, FWAction
import pathlib, logging
import pycp2k, cp2kparser
import ase, ase.io
from critcatworks.database import atoms_dict_to_ase, ase_to_atoms_dict
from critcatworks.database.extdb import update_simulations_collection
import copy

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


        simulation = fw_spec["simulation"]
        atoms_dict = simulation["atoms"]
        atoms = atoms_dict_to_ase(atoms_dict)
        cell_size = atoms.get_cell()
        


        calc.CP2K_INPUT.FORCE_EVAL_list[0].SUBSYS.CELL.A = "[angstrom] " + str(cell_size[0][0]) + " " + str(cell_size[0][1]) + " " + str(cell_size[0][2]) 
        calc.CP2K_INPUT.FORCE_EVAL_list[0].SUBSYS.CELL.B = "[angstrom] " + str(cell_size[1][0]) + " " + str(cell_size[1][1]) + " " + str(cell_size[1][2]) 
        calc.CP2K_INPUT.FORCE_EVAL_list[0].SUBSYS.CELL.C = "[angstrom] " + str(cell_size[2][0]) + " " + str(cell_size[2][1]) + " " + str(cell_size[2][2]) 
        calc.CP2K_INPUT.FORCE_EVAL_list[0].SUBSYS.TOPOLOGY.Coord_file_name = "structure.xyz"
        calc.working_directory = str(target_path)
        logging.debug("working_directory: " + str(calc.working_directory))
        calc.project_name = "gopt"
        calc.write_input_file()
        logging.info("cp2k input file written TO" + calc.project_name + ".inp")


        input_string = calc.get_input_string()
        update_spec = fw_spec
        update_spec["simulation"]["inp"]["input_string"] = input_string
        update_spec["simulation"]["inp"]["path"] = str(target_path)
        #pass_spec = fw_spec
        #print("dummy outputs generated")
        #fw_spec.pop("_category")
        #fw_spec.pop("name")
        #detours = Firework([CP2KRunTask(target_path=target_path, calc_id = calc_id, n_max_restarts = n_max_restarts)], 
        #    spec = {'_category' : "dft", 'name' : 'CP2KRunTask', "n_restarts" : 0},
        #    name = 'CP2KRunWork')
        #return FWAction(update_spec = fw_spec, detours = detours)
        #return FWAction(update_spec = fw_spec)
        return FWAction()

@explicit_serialize
class CP2KRunTask(FiretaskBase):
    """ 
    Task to run CP2K calculations.

    Args:
        None
    """

    _fw_name = 'CP2KRunTask'
    required_params = ['target_path', 'calc_id', "n_max_restarts", 'skip_dft']
    optional_params = []

    def run_task(self, fw_spec):
        target_path = self["target_path"]
        calc_id = self["calc_id"]
        n_max_restarts = self["n_max_restarts"]
        skip_dft = self["skip_dft"]
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
            if skip_dft:
                logging.warning("DFT calculation is skipped. Switch skip_dft to False!")
            else:
                subprocess.call(command_list, shell = False)
                print("run done")

        fw_spec.pop("_category")
        fw_spec.pop("name")
        # remove detour 
        #detours = Firework([CP2KAnalysisTask(target_path=target_path, calc_id = calc_id, n_max_restarts = n_max_restarts, skip_dft = skip_dft)], 
        #    spec = {'_category' : "lightweight", 'name' : 'CP2KAnalysisTask', 
        #        "n_restarts" : n_restarts, "simulation" : fw_spec["simulation"],
        #        "extdb_connect" : fw_spec["extdb_connect"]},
        #    name = 'CP2KAnalysisWork')
        return FWAction(update_spec = fw_spec) #, detours = detours)


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
        skip_dft = self["skip_dft"]
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
        if output_state == "no_output": # or output_state == "incorrect_termination":
            pass
        else:
            nwarnings = preparse_results.get('nwarnings', 0)
            is_walltime_exceeded = preparse_results.get('exceeded_walltime', False)

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
            cell = results["simulation_cell"][-1]

            input_filename = results["x_cp2k_input_filename"]
            with open(target_path + "/" + input_filename, "r") as f:
                input_string = f.read()

            if is_converged:
                total_energy = frame_sequence_potential_energy[-1]
            else:
                output_state = "not_converged"
                logging.info("CP2K not converged")
                print("CP2K not converged")
                # not converged energy is stored
                total_energy = frame_sequence_potential_energy[-1]

        # restart
        if output_state == "no_output" or output_state == "incorrect_termination" or output_state == "not_converged":
            if fw_spec["n_restarts"] < n_max_restarts:
                fw_spec["n_restarts"] += 1
                detours =  rerun_cp2k(target_path, calc_id, n_max_restarts, n_restarts = int(n_restarts) + 1, 
                    simulation = fw_spec["simulation"], skip_dft = skip_dft, extdb_connect = fw_spec["extdb_connect"])

                fw_spec.pop("_category")
                fw_spec.pop("name")
                return FWAction(update_spec = fw_spec, detours = detours)
            else:
                is_defused = True   
        else:
            is_defused = False 

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
    
            # conversion factor cp2kparser uses other units! WARNING, might change in the future! conversion always back to Angstrom.
            relaxed_structure = ase.Atoms(symbols = atom_labels[-1], 
                positions = atom_positions[-1] * 10e9, cell = cell * 10e9)
    
            atoms_dict = ase_to_atoms_dict(relaxed_structure)
    
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
            atoms_dict = fw_spec["simulation"]["atoms"]
            result_dict = {"is_converged" : 0, "output_state" : output_state}
            input_string = "no output or incorrect termination of CP2K calculation"

        ##
        # get source simulation
        source_simulation = fw_spec["simulation"]

        # update external database
        dct = copy.deepcopy(source_simulation)
        # calculation originated from this:
        dct["source_id"] = calc_id
        dct["atoms"] = atoms_dict ### !!!
        dct["operations"] = ["cp2k"]
        dct["output"] = result_dict # might still be missing some output
        dct["inp"]["path"] = str(target_path)
        dct["inp"]["input_string"] = str(input_string)
        logging.info("simulation after analysis")

        simulation = update_simulations_collection(extdb_connect = fw_spec["extdb_connect"], **dct)

        logging.info(simulation)

        # update internal workflow data
        simulation_id = simulation["_id"]

        # update temp workflow data
        mod_spec =[
            #{"_set" : {"simulations->" + str(simulation_id) : simulation }},
            {"_push" : {"temp->" + "analysis_ids" : simulation_id }},
            ]

        logging.info("pushing new ids to analysis_ids")

        fw_spec.pop("_category")
        fw_spec.pop("name")
        # update spec no longer necessary
        #return FWAction(update_spec = fw_spec, mod_spec=mod_spec)
        return FWAction(mod_spec=mod_spec, defuse_children = is_defused)


def setup_cp2k(template, target_path, calc_id, simulation, name = "cp2k_run_id", n_max_restarts = 4,
        skip_dft = False, extdb_connect = {}):
    setup_task = CP2KSetupTask(template = template,
                    target_path = target_path,
                    calc_id = calc_id,
                    name = name,
                    n_max_restarts = n_max_restarts,
                    )

    run_task = CP2KRunTask(target_path = target_path, 
        calc_id = calc_id, 
        n_max_restarts = n_max_restarts,
        skip_dft = skip_dft)
    
    # add analysis as an additional Firework as a child. Ensures
    # that Firework exists even if CP2K fails
    analysis_task = CP2KAnalysisTask(target_path = target_path, 
        calc_id = calc_id, 
        n_max_restarts = n_max_restarts,
        skip_dft = skip_dft)

    
    fw1 = Firework([setup_task,run_task], 
        spec = {'_category' : "dft", 'name' : 'CP2KWork', 
            "n_restarts" : 0, "simulation" : simulation,
            "extdb_connect" : extdb_connect,
            "_priority": 10,},
        name = 'CP2KWork')

    fw2 = Firework([analysis_task], 
        spec = {'_category' : "lightweight", 'name' : 'CP2KAnalysisTask', 
            "n_restarts" : 0,
            "simulation" : simulation,
            "extdb_connect" : extdb_connect,
            "_allow_fizzled_parents" : True,
            "_priority": 8,},
        name = 'CP2KAnalysisWork',
        parents = [fw1]) 
    #return Workflow([fw1,fw2], {fw1 : [fw2]})
    #return [fw1,fw2], {fw1 : [fw2]}
    return [fw1,fw2], {fw1 : [fw2]}


def rerun_cp2k(target_path, calc_id, n_max_restarts, n_restarts, 
        simulation, skip_dft = False, extdb_connect = {}):    
    """
    Run and analyse
    """
    fw1 = Firework([CP2KRunTask(target_path=target_path, calc_id = calc_id, 
        n_max_restarts = n_max_restarts, skip_dft = skip_dft)], 
        spec = {'_category' : "dft", 'name' : 'CP2KRerunWork', 
            'n_restarts' : int(n_restarts), "simulation" : simulation,
            "_priority": 10,
            "extdb_connect" : extdb_connect },
        name = 'CP2KRerunWork')

    # add analysis as an additional Firework as a child. Ensures
    # that Firework exists even if CP2K fails
    analysis_task = CP2KAnalysisTask(target_path = target_path, 
        calc_id = calc_id, 
        n_max_restarts = n_max_restarts,
        skip_dft = skip_dft)

    fw2 = Firework([analysis_task], 
        spec = {'_category' : "lightweight", 'name' : 'CP2KAnalysisTask', 
            'n_restarts' : int(n_restarts), "simulation" : simulation,
            "extdb_connect" : extdb_connect,
            "_priority": 8,
            "_allow_fizzled_parents" : True},
        name = 'CP2KAnalysisWork',
        parents = [fw1]) 
    return Workflow([fw1,fw2],)
    #return Workflow([fw1,fw2], {fw1 : [fw2]})
    #return [fw1,fw2], {fw1 : [fw2]}

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
