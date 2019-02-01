# Setup of folders
# ChunkCalculations
from fireworks import Firework, FWorker, LaunchPad, Workflow
import os,time, re, glob, sys, subprocess
from fireworks import explicit_serialize, FiretaskBase, FWAction
import pathlib, logging
import ase, ase.io
from critcatworks.database import atoms_dict_to_ase
from critcatworks.dft import setup_cp2k

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

            ase.io.write(structure_folder_path + "/" + "structure.xyz", atoms)
        
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
    required_params = ['template_path', 'target_path',  'name', 'n_max_restarts']
    optional_params = ['chunk_size', 'simulation_method']

    def run_task(self, fw_spec):
        template_path = self["template_path"]
        target_path = self["target_path"]
        chunk_size = self.get("chunk_size", -1)
        simulation_method = self.get("simulation_method", "cp2k")
        name = self["name"]
        n_max_restarts = self["n_max_restarts"]
        calc_paths = fw_spec["calc_paths"]
        #fps_ranking = fw_spec["fps_ranking"]
        calc_ids = fw_spec["calc_ids"]


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



        #ranked_ids = fps_ranking[n_calcs_started : n_calcs_started+chunk_size]

        detours = []
        for calc_id in calc_ids:
            logging.info(calc_id)
            target_path = calc_paths[int(calc_id)]

            if simulation_method == "cp2k":
                # create detour to setup cp2k calculation
                new_fw = setup_cp2k(template_path = template_path,
                    target_path = target_path,
                    calc_id = calc_id,
                    name = name,
                    n_max_restarts = n_max_restarts)
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


def chunk_calculations(template_path, target_path, chunk_size = None, name = "cp2k_run_id", n_max_restarts = 4, simulation_method = "cp2k"):
    firetask1  = ChunkCalculationsTask(
        template_path = template_path,
        target_path = target_path,
        chunk_size = chunk_size,
        name = name,
        n_max_restarts = n_max_restarts,
        simulation_method = simulation_method,
        )
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'ChunkCalculationsTask'},
                     name = 'ChunkCalculationsWork')
    return fw