from fireworks import LaunchPad, Workflow

import pathlib
import os,time

# internal modules
from critcatworks.database.update import initialize_workflow_data
from critcatworks.database.format import ase_to_atoms_dict
from critcatworks.database import start_from_structures, start_from_database, read_structures
from critcatworks.structure import get_monodentate_unique_adsites
from critcatworks.dft import setup_folders, chunk_calculations

from critcatworks.structure import update_converged_data
from critcatworks.ml import get_mae, check_convergence


def get_uniquemolsites_workflow(template_path, username, password, 
        worker_target_path = None, structures = None, extdb_ids = None,
        source_path  = None, reference_energy=0.0, adsorbate = {}, 
        adsite_types = ["top", "bridge", "hollow"], threshold = 0.001, n_max_restarts = 1,
        skip_dft = False, extdb_connect = {}):
    """ 
    Workflow to determine the adsorption sites and energies of a set of
    nanocluster structures. The adsorption sites are determined by the 
    python package cluskit and then ranked by farthest point sampling
    based on their structural local dissimilarity.
    Only sites which are more dissimilar than the given threshold are 
    computed.
    The adsorption energy is determined by a simulation code (e.g. CP2K).

    Args:
        template_path (str) :   absolute path to input file for calculations. 
                                It works as a template which is later modified by the
                                simulation-specific Firework.
        username (str) :        user who executed the workflow
        password (str) :        password for user to upload to the database
        worker_target_path (str) :  absolute path on computing resource 
                                    directory needs to exist
        structures (list) : list of ase.Atoms objects from where the workflow is started.
        extdb_ids (list) :  unique identifiers of the simulations collection which
                            are used to start the workflow
        source_path (str) : absolute path on the computing resource to the directory 
                            where to read the structures from        
        reference_energy (float) :  reference energy for the adsorbate. Can be the
                                    total energy of the isolated adsorbate molecule
                                    or a different reference point
        adsorbate (dict) :  adsorbed molecule as atoms dict. Contains an "X" dummy atom
                            which indicates the anchor point to the nanocluster
        adsite_types (list) :   adsorption site types, can contain any combination of
                                "top", "bridge", "hollow"
        threshold (float) :     threshold of similarity metric between the local structures
                                of the adsorption sites. Only sites which are more dissimilar 
                                than the given threshold are computed
        n_max_restarts (int)  : number of times the calculation is restarted upon failure
        skip_dft (bool) :   If set to true, the simulation step is skipped in all
                            following simulation runs. Instead the structure is returned unchanged.
        extdb_connect (dict):   dictionary containing the keys host,
                                username, password, authsource and db_name.
        
    Returns:
        fireworks.Workflow : molsinglesites Fireworks Workflow object
    """

    with open (template_path, "r") as f:
        template = f.read()

    #translate adsorbate molecule to json format
    adsorbate_dict = ase_to_atoms_dict(adsorbate)

    #FireWork: Initialize workflow with workflow_id from external database
    parameters = {
        "template" : template,
        "template_path" : template_path,
        "worker_target_path" : worker_target_path,
        "extdb_ids" : extdb_ids,
        "source_path" : source_path,
        "reference_energy" : reference_energy,
        "adsorbate" : adsorbate_dict,
        "adsite_types" : adsite_types,
        "descriptor" : "soap",
        "descriptor_params" : {"nmax" : 9, "lmax" :6, "rcut" : 5.0, 
            "crossover" : True, "sparse" : False},
        "simulation_method" : "cp2k",
        "threshold" : threshold,
        "n_max_restarts" : n_max_restarts,
        "workflow_type" : "molsinglesites",
        }

    fw_init = initialize_workflow_data(username, password, parameters, 
        name = "UNNAMED", workflow_type = "uniquemolsites",
        extdb_connect = extdb_connect)

    # FireWork: Read nanocluster structures and initialise a database
    # object containing set information
    if structures != None:
        jsonified_structures = []
        for atoms in structures:
            atoms_dict = ase_to_atoms_dict(atoms)
            jsonified_structures.append(atoms_dict) 
        fw_get_structures = start_from_structures(jsonified_structures)
    elif extdb_ids != None:
        fw_get_structures = start_from_database(extdb_ids)
    elif source_path != None:
        fw_get_structures = read_structures(source_path)
    else:
        raise ValueError('structures, extdb_ids or source_path contain no entries!')

    # FireWork: Determine adsites and add to database
    fw_get_adsites = get_monodentate_unique_adsites(
        reference_energy= reference_energy, 
        adsorbate = adsorbate_dict, 
        adsite_types = adsite_types,
        descriptor = "soap",
        descriptor_params = {"nmax" : 9, "lmax" :6, "rcut" : 5.0, 
            "crossover" : True, "sparse" : False},
        threshold = threshold,
        )

    # Firework: setup folders for DFT calculations
    fw_setup_folders = setup_folders(target_path = worker_target_path, name = "cp2k_uniquemolsites_id")

    # add above Fireworks with links
    workflow_list = [fw_init,
        fw_get_structures, 
        fw_get_adsites, 
        fw_setup_folders,
        ]

    links_dict = {
            fw_init: [fw_get_structures],
            fw_get_structures: [fw_get_adsites], 
            fw_get_adsites : [fw_setup_folders],
            }

    # FireWork: setup, run and extract DFT calculation
    # (involves checking for errors in DFT and rerunning)
    fw_chunk_calculations = chunk_calculations(template = template, target_path = worker_target_path, 
        n_max_restarts = n_max_restarts, simulation_method = "cp2k",
        skip_dft = skip_dft)
    workflow_list.append(fw_chunk_calculations)
    links_dict[fw_setup_folders] = [fw_chunk_calculations]

    # FireWork: update database, 
    # (includes reading relaxed structure and energy)
    fw_update_converged_data = update_converged_data(chunk_size=-1)
    workflow_list.append(fw_update_converged_data)
    links_dict[fw_chunk_calculations] =[fw_update_converged_data]

    wf = Workflow(workflow_list, links_dict, name = "uniquemolsites")
    return wf
