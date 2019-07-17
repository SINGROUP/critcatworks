# Setup of folders
# ChunkCalculations
from fireworks import Firework, FWorker, LaunchPad, Workflow
import os,time, re, glob, sys, subprocess
from fireworks import explicit_serialize, FiretaskBase, FWAction
import pathlib, logging
import ase, ase.io
from critcatworks.database import atoms_dict_to_ase, ase_to_atoms_dict
from critcatworks.dft.cp2k import setup_cp2k
from critcatworks.database.extdb import fetch_simulations

@explicit_serialize
class StructureFolderTask(FiretaskBase):
    """ 
    Task to setup folders with xyz structures.

    Args:
        target_path (str) : absolute path to the target directory 
                            (needs to exist) on the computing resource.
        name (str) :        individual calculation folder name 
                            is prefixed with the given string

    Returns:
        FWAction : Firework action, updates fw_spec
    """
    _fw_name = 'StructureFolderTask'
    required_params = ['target_path', 'name']
    optional_params = []

    def run_task(self, fw_spec):
        target_path = self['target_path']
        prefix = self['name']
        time_str = time.strftime("%Y-%m-%d-%H-%M")

        parent_folder_name = 'cp2k_calculations_' + time_str
        parent_folder_path = target_path + "/" + parent_folder_name
        calc_ids = fw_spec["temp"]["calc_ids"]
        simulations = fetch_simulations(fw_spec["extdb_connect"], calc_ids) 

        if not os.path.exists(parent_folder_path):
            os.makedirs(parent_folder_path)

        calc_paths = []
        #iterating over available structures
        for idx, calc_id in enumerate(calc_ids):

            atoms_dict = simulations[str(calc_id)]["atoms"]
            atoms = atoms_dict_to_ase(atoms_dict)
            structure_folder = prefix + '_' + str(idx)

            structure_folder_path = parent_folder_path + "/" + structure_folder
            calc_paths.append(structure_folder_path)

            if not os.path.exists(structure_folder_path):
                os.makedirs(structure_folder_path)

            ase.io.write(structure_folder_path + "/" + "structure.xyz", atoms)
        
        update_spec = fw_spec


        update_spec["temp"]["calc_paths"] = calc_paths
        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec = update_spec)


@explicit_serialize
class ChunkCalculationsTask(FiretaskBase):
    """ 
    Create Fireworks with new calculations to setup and run

    Args:
        template (str)    : input file for calculations represented as string. 
                            It works as a template which is later modified by the
                            simulation-specific Firework.
        target_path (str) : absolute path to the target directory 
                            (needs to exist) on the computing resource.
        name (str) :        individual calculation folder name 
                            is prefixed with the given string
        n_max_restarts (int)  : number of times the calculation is restarted upon failure
        chunk_size (int) :  number of calculations to be run simulataneously. Default -1
                            means all calculations are run at once.
        simulation_method (str) :   Specifies which simulation code to use.
                                    Currently, only CP2K is implemented.
        skip_dft (bool) :   If set to true, the simulation step is skipped in all
                            following simulation runs. Instead the structure is returned unchanged.
        is_safeguard (bool) : if False, the workflow is not paused when not all CP2K jobs
                               converge properly after the maximum number of restarts.
    Returns:
        FWAction :  Firework action, updates fw_spec, 
                    creates new Fireworks as detours from workflow
    """

    _fw_name = 'ChunkCalculationsTask'
    required_params = ['template', 'target_path',  'name', 'n_max_restarts']
    optional_params = ['chunk_size', 'simulation_method', 'skip_dft', 'is_safeguard']

    def run_task(self, fw_spec):
        template = self["template"]
        target_path = self["target_path"]
        chunk_size = self.get("chunk_size", -1)
        is_safeguard = self.get("is_safeguard", True)
        simulation_method = self.get("simulation_method", "cp2k")
        name = self["name"]
        n_max_restarts = self["n_max_restarts"]
        skip_dft = self["skip_dft"]
        calc_paths = fw_spec["temp"]["calc_paths"]
        calc_ids = fw_spec["temp"]["calc_ids"]
        simulations = fetch_simulations(fw_spec["extdb_connect"], calc_ids)

        # define what chunk to run
        try:
            n_calcs_started = fw_spec["n_calcs_started"]
        except KeyError:
            logging.info("Starting Chunk")
            n_calcs_started = 0

        if chunk_size == -1:
            pass
        else:
            calc_ids = calc_ids[n_calcs_started : n_calcs_started+chunk_size]
            calc_paths = calc_paths[n_calcs_started : n_calcs_started+chunk_size]

        wfs = []
        links_dict = {}
        for idx, calc_id in enumerate(calc_ids):
            logging.info(calc_id)
            target_path = calc_paths[int(idx)]

            if simulation_method == "cp2k":
                # create detour to setup cp2k calculation
                simulation = simulations[str(calc_id)]
                new_fw, links = setup_cp2k(template = template,
                    target_path = target_path,
                    calc_id = calc_id,
                    name = name,
                    n_max_restarts = n_max_restarts,
                    simulation = simulation,
                    skip_dft = skip_dft,
                    is_safeguard = is_safeguard,
                    extdb_connect = fw_spec["extdb_connect"]
                    )
                wfs.extend(new_fw)
                links_dict.update(links_dict)

            else:
                logging.warning("WARNING! " + str(simulation_method) + " unknown simulation method.")
        detours = Workflow(wfs, links_dict)

        update_spec = fw_spec
        update_spec["n_calcs_started"] = n_calcs_started + chunk_size
        update_spec.pop("_category")
        update_spec.pop("name")

        return FWAction(update_spec = update_spec, detours = detours)


def setup_folders(target_path, name = "cp2k_run_id",):
    """ 
    Creates folders with xyz structures.

    Args:
        target_path (str) : absolute path to the target directory 
                            (needs to exist) on the computing resource.
        name (str) :        individual calculation folder name 
                            is prefixed with the given string

    Returns:
        Firework : StructureFolderWork Firework
    """
    firetask1  = StructureFolderTask(
        target_path = target_path,
        name = name)
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'StructureFolderTask'},
             name = 'StructureFolderWork')
    return fw


def chunk_calculations(template, target_path, chunk_size = -1, name = "cp2k_run_id", n_max_restarts = 4, simulation_method = "cp2k", 
    skip_dft = False, is_safeguard = True):
    """ 
    Create Fireworks with new calculations to setup and run.

    Args:
        template (str)    : input file for calculations represented as string. 
                            It works as a template which is later modified by the
                            simulation-specific Firework.
        target_path (str) : absolute path to the target directory 
                            (needs to exist) on the computing resource.
        name (str) :        individual calculation folder name 
                            is prefixed with the given string
        n_max_restarts (int)  : number of times the calculation is restarted upon failure
        chunk_size (int) : det :    number of calculations to be run simulataneously. Default -1
                                    means all calculations are run at once.
        simulation_method (str) :   Specifies which simulation code to use.
                                    Currently, only CP2K is implemented.
        skip_dft (bool) :   If set to true, the simulation step is skipped in all
                            following simulation runs. Instead the structure is returned unchanged.
        is_safeguard (bool) : if False, the workflow is not paused when not all CP2K jobs
                               converge properly after the maximum number of restarts.

    Returns:
                Firework : StructureFolderWork Firework,
                           creates new Fireworks as detours from workflow
    """
    firetask1  = ChunkCalculationsTask(
        template = template,
        target_path = target_path,
        chunk_size = chunk_size,
        name = name,
        n_max_restarts = n_max_restarts,
        simulation_method = simulation_method,
        skip_dft = skip_dft,
        is_safeguard = is_safeguard,
        )
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'ChunkCalculationsTask'},
                     name = 'ChunkCalculationsWork')
    return fw
