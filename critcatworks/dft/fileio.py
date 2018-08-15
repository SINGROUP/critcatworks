from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time
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
    required_params = ['template', 'target_path', 'chunk_size', 'name']
    optional_params = []

    def run_task(self, fw_spec):
        template = self["template"]
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
            new_fw = Firework([CP2KSetupTask(template = template,
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
    required_params = ['template', 'target_path', 'ranked_id']
    optional_params = ['name']

    def run_task(self, fw_spec):
        
        logging.info("CP2KSetupTask not implemented yet")
        
        prefix = self.get("name", "cp2k_run_id")      
        template_path = self["template"]
        target_path = self["target_path"]
        ranked_id = self["ranked_id"]

        # read template

        inpparser = pycp2k.CP2KInputParser()
        calc = inpparser.parse(template_path)

        logging.info("info about input parser")
        logging.info(inpparser)

        logging.info("cp2k info storage \n")
        logging.info(inpparser.storage_obj)


        logging.debug("target_path", target_path)

        calc.working_directory = str(target_path)
        logging.debug("working_directory", calc.working_directory)
        calc.project_name = "gopt"
        calc.write_input_file()

        logging.info("cp2k input file written TO" + calc.project_name + ".inp")
        #pass_spec = fw_spec
        #print("dummy outputs generated")

        detours = Firework([CP2KRunTask(target_path=target_path, ranked_id = ranked_id)])
        return FWAction(detours = detours)



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
        logging.info("Running CP2K not implemented yet")
        detours = Firework([CP2KAnalysisTask(target_path=target_path, ranked_id = ranked_id)])
        return FWAction(detours = detours)


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
        total_energy = ranked_id
        logging.info(ranked_id)

        atoms_dict = "NotYetImplemented"

        mod_spec =[
            {'_set' : {'adsorbate_energies_dict->' + str(ranked_id): total_energy}},
            {'_set' : {'relaxed_structure_dict->' + str(ranked_id): atoms_dict}},
            ]
        return FWAction(mod_spec=mod_spec)



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
       logging.info("DFTSetupTask not implemented yet")
       logging.info("forwarding information")
       pass_spec = fw_spec
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
    logging.info(atoms.get_potential_energy())
    return
###############################################################################################