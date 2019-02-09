# Setup of folders
# ChunkCalculations
from fireworks import Firework, FWorker, LaunchPad, Workflow
import os,time, re, glob, sys, subprocess
from fireworks import explicit_serialize, FiretaskBase, FWAction
import pathlib, logging
import ase, ase.io
from critcatworks.database import atoms_dict_to_ase, ase_to_atoms_dict
from critcatworks.dft.cp2k import setup_cp2k

@explicit_serialize
class StructureFolderTask(FiretaskBase):
    """ 
    Task to setup folders with xyz structures.

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
        simulations = fw_spec["simulations"]
        calc_ids = fw_spec["temp"]["calc_ids"]

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
        None
    """

    _fw_name = 'ChunkCalculationsTask'
    required_params = ['template', 'target_path',  'name', 'n_max_restarts']
    optional_params = ['chunk_size', 'simulation_method']

    def run_task(self, fw_spec):
        template = self["template"]
        target_path = self["target_path"]
        chunk_size = self.get("chunk_size", -1)
        simulation_method = self.get("simulation_method", "cp2k")
        name = self["name"]
        n_max_restarts = self["n_max_restarts"]
        skip_dft = self["skip_dft"]
        calc_paths = fw_spec["temp"]["calc_paths"]
        calc_ids = fw_spec["temp"]["calc_ids"]


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

        detours = []
        for idx, calc_id in enumerate(calc_ids):
            logging.info(calc_id)
            target_path = calc_paths[int(idx)]

            if simulation_method == "cp2k":
                # create detour to setup cp2k calculation
                simulation = fw_spec["simulations"][str(calc_id)]
                new_fw = setup_cp2k(template = template,
                    target_path = target_path,
                    calc_id = calc_id,
                    name = name,
                    n_max_restarts = n_max_restarts,
                    simulation = simulation,
                    skip_dft = skip_dft,
                    extdb_connect = fw_spec["extdb_connect"]
                    )
                detours.append(new_fw)

            else:
                logging.warning("WARNING! " + str(simulation_method) + " unknown simulation method.")


        update_spec = fw_spec
        update_spec["n_calcs_started"] = n_calcs_started + chunk_size
        update_spec.pop("_category")
        update_spec.pop("name")

        return FWAction(update_spec = update_spec, detours = detours)


def setup_folders(target_path, name = "cp2k_run_id",):
    firetask1  = StructureFolderTask(
        target_path = target_path,
        name = name)
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'StructureFolderTask'},
             name = 'StructureFolderWork')
    return fw


def chunk_calculations(template, target_path, chunk_size = -1, name = "cp2k_run_id", n_max_restarts = 4, simulation_method = "cp2k", skip_dft = False):
    firetask1  = ChunkCalculationsTask(
        template = template,
        target_path = target_path,
        chunk_size = chunk_size,
        name = name,
        n_max_restarts = n_max_restarts,
        simulation_method = simulation_method,
        skip_dft = skip_dft,
        )
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'ChunkCalculationsTask'},
                     name = 'ChunkCalculationsWork')
    return fw
