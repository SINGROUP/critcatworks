from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time, re, glob
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import pathlib, logging
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
        return FWAction(update_spec = update_spec)

@explicit_serialize
class ChunkCalculationsTask(FiretaskBase):
    """ 
    Create Fireworks with new calculations to setup and run

    Args:
        None
    """

    _fw_name = 'ChunkCalculationsTask'
    required_params = ['template_path', 'target_path', 'chunk_size', 'name']
    optional_params = []

    def run_task(self, fw_spec):
        template_path = self["template_path"]
        target_path = self["target_path"]
        chunk_size = self["chunk_size"]
        name = self["name"]

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
                name = name)])
            detours.append(new_fw)


        update_spec = fw_spec
        update_spec["n_calcs_started"] = n_calcs_started + chunk_size

        return FWAction(update_spec = update_spec, detours = detours)


@explicit_serialize
class CP2KSetupTask(FiretaskBase):
    """ 
    Task to setup DFT calculations.

    Args:
        None
    """

    _fw_name = 'CP2KSetupTask'
    required_params = ['template_path', 'target_path', 'ranked_id']
    optional_params = ['name']

    def run_task(self, fw_spec):
        
        logging.info("CP2KSetupTask not implemented yet")
        
        prefix = self.get("name", "cp2k_run_id")      
        template_path = self["template_path"]
        target_path = self["target_path"]
        ranked_id = self["ranked_id"]

        # read template
        cp2kinput = glob.glob(template_path + "/" + "*inp")[0]

        inpparser = pycp2k.CP2KInputParser()
        calc = inpparser.parse(cp2kinput)

        logging.info("info about input parser")
        logging.info(inpparser)

        logging.info("cp2k info storage \n")
        logging.info(inpparser.storage_obj)


        logging.debug("target_path", target_path)
        
        calc.CP2K_INPUT.FORCE_EVAL_list[0].SUBSYS.CELL.Abc = "[angstrom] 20 20 20"

        calc.working_directory = str(target_path)
        logging.debug("working_directory", calc.working_directory)
        calc.project_name = "gopt"
        calc.write_input_file()

        logging.info("cp2k input file written TO" + calc.project_name + ".inp")
        #pass_spec = fw_spec
        #print("dummy outputs generated")

        detours = Firework([CP2KRunTask(target_path=target_path, ranked_id = ranked_id)])
        return FWAction(update_spec = fw_spec, detours = detours)


@explicit_serialize
class CP2KRunTask(FiretaskBase):
    """ 
    Task to run CP2K calculations.

    Args:
        None
    """

    _fw_name = 'CP2KRunTask'
    required_params = ['target_path', 'ranked_id']
    optional_params = []

    def run_task(self, fw_spec):
        target_path = self["target_path"]
        ranked_id = self["ranked_id"]
        logging.info("Running CP2K not implemented yet. Creating dummy outputs")

        with open(target_path + "/fake.out", "w") as f:
            f.write("WARNING")
            f.write("total_energy:        42.42")

        detours = Firework([CP2KAnalysisTask(target_path=target_path, ranked_id = ranked_id)])
        return FWAction(update_spec = fw_spec, detours = detours)


@explicit_serialize
class CP2KAnalysisTask(FiretaskBase):
    """ 
    Task to analyse CP2K calculations.

    Args:
        None
    """

    _fw_name = 'CP2KAnalysisTask'
    required_params = ['target_path', 'ranked_id']
    optional_params = []

    def run_task(self, fw_spec):
        target_path = self["target_path"]
        ranked_id = self["ranked_id"]
        
        logging.info("Analysis of CP2K output not implemented yet")


        # read output
        cp2koutput = glob.glob(target_path + "/" + "*out")[0]
        total_energy = 42.42 + float(ranked_id / 100.0)
        with open(cp2koutput) as origin_file:
            for line in origin_file:
                print(line)
                is_toten = re.findall(r'total_energy:', line)
                if is_toten:
                    total_energy = line.split(None)[1]
                    print(total_energy)

                is_warning = re.findall(r'Warning:', line, flags=re.IGNORECASE)

                is_converged = re.findall(r'Warning:', line, flags=re.IGNORECASE)

                is_ended = re.findall(r'PROGRAM ENDED AT', line)


        logging.info(ranked_id)

        atoms_dict = "NotYetImplemented"

        mod_spec =[
            {'_set' : {'adsorbate_energies_dict->' + str(ranked_id): float(total_energy)}},
            {'_set' : {'relaxed_structure_dict->' + str(ranked_id): atoms_dict}},
            ]
        return FWAction(update_spec = fw_spec, mod_spec=mod_spec)


def setup_cp2k(template_path, target_path, chunk_size, name = "cp2k_run_id",):
    firetask1  = ChunkCalculationsTask(
        template_path = template_path,
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
