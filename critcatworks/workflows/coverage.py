from fireworks import LaunchPad, Workflow

import pathlib
import os,time

# internal modules
from critcatworks.database.update import initialize_workflow_data
from critcatworks.database.format import ase_to_atoms_dict
from critcatworks.database import start_from_structures, start_from_database, read_structures
from critcatworks.structure import update_converged_data
from critcatworks.dft import setup_folders, chunk_calculations
from critcatworks.structure import get_per_type_coverage, eliminate_pairs, eliminate_closest


def get_coverage_workflow(template_path, username, password, 
        worker_target_path = None, structures = None, extdb_ids = None,
        source_path  = None, reference_energy=0.0,
        adsorbate_name='H',max_iterations = 10000,  
        adsite_types = ["top", "bridge", "hollow"], n_max_restarts = 1, 
        skip_dft = False, is_safeguard = True, bond_length = 1.4, n_remaining = "",
        extdb_connect = {}):
    """
    Workflow to determine a stable coverage of a nanocluster with single adsorbate atoms. As a first step, 
    adsorbates are put on top, bridge and hollow sites. Once the structure is relaxed by DFT,
    formed adsorbate molecules (pairs of atoms) are replaced by a single adsorbate.
    The procedure is repeated until no adsorbate molecules form.
    
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
        adsorbate_name (str) : element symbol of the adsorbed atom
        max_iterations (int) : maximum number of iterations in the workflow
        adsite_types (list) :   adsorption site types, can contain any combination of
                                "top", "bridge", "hollow"
        n_max_restarts (int)  : number of times the calculation is restarted upon failure
        skip_dft (bool) :   If set to true, the simulation step is skipped in all
                            following simulation runs. Instead the structure is returned unchanged.
        is_safeguard (bool) : if False, the workflow is not paused when not all simulation jobs
                               converge properly after the maximum number of restarts.
        bond_length (float) :   distance in angstrom under which two adsorbed atoms are 
                                considered bound, hence too close
        n_remaining (int) : number of adsorbates which should remain after the
                            first pre-DFT pruning of the adsorbate coverage
        extdb_connect (dict):   dictionary containing the keys host,
                                username, password, authsource and db_name.
        
    Returns:
        fireworks.Workflow : coverage Fireworks Workflow object
    """
    
    with open (template_path, "r") as f:
        template = f.read()

    #FireWork: Initialize workflow with workflow_id from external database
    parameters = {
        "template" : template,
        "template_path" : template_path,
        "worker_target_path" : worker_target_path,
        "extdb_ids" : extdb_ids,
        "source_path" : source_path,
        "reference_energy" : reference_energy,
        "max_iterations" : max_iterations,
        "adsorbate_name" : adsorbate_name,
        "adsite_types" : adsite_types,
        "descriptor" : "soap",
        "descriptor_params" : {"nmax" : 9, "lmax" :6, "rcut" : 5.0, 
            "crossover" : True, "sparse" : False},
        "simulation_method" : "cp2k",
        "n_max_restarts" : n_max_restarts,
        "workflow_type" : "pertype_coverage",
        }

    fw_init = initialize_workflow_data(username, password, parameters, 
        name = "UNNAMED", workflow_type = "coverage",
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
    # create structure with coverage
    fw_get_per_type_coverage = get_per_type_coverage(
        reference_energy=reference_energy, 
        adsorbate_name='H', 
        adsite_types = ["top", "bridge", "hollow"],
        descriptor = "soap",
        descriptor_params = {"nmax" : 9, "lmax" :6, "rcut" : 5.0, 
            "crossover" : True, "sparse" : False},
        )

    # FireWork: before running DFT eliminate too close adsorbates 
    # eliminate adsorbate pairs too close
    if n_remaining:
        fw_eliminate_pairs = eliminate_closest(adsorbate_name = adsorbate_name, n_remaining = n_remaining)
    else:
        fw_eliminate_pairs = eliminate_pairs(adsorbate_name = adsorbate_name, bond_length = bond_length)

    # add above Fireworks with links
    workflow_list = [fw_init,
        fw_get_structures, 
        fw_get_per_type_coverage,
        fw_eliminate_pairs,
        ]

    links_dict = {
            fw_init : [fw_get_structures],
            fw_get_structures: [fw_get_per_type_coverage], 
            fw_get_per_type_coverage : [fw_eliminate_pairs],
            }

    ### loop starts ###
    for i in range(max_iterations):
        # Firework: setup folders for DFT calculations,
        fw_setup_folders = setup_folders(target_path = worker_target_path, name = "cp2k_coverage_iter_" + str(i))
        workflow_list.append(fw_setup_folders)
        links_dict[fw_eliminate_pairs] = [fw_setup_folders]

        # FireWork: setup, run and extract DFT calculation
        # (involves checking for errors in DFT and rerunning)
        fw_chunk_calculations = chunk_calculations(template = template, 
            target_path = worker_target_path, 
            chunk_size = -1, n_max_restarts = n_max_restarts, 
            simulation_method = "cp2k",
            name = "cp2k_coverage_iter_" + str(i), 
            skip_dft = skip_dft, is_safeguard = is_safeguard)
        workflow_list.append(fw_chunk_calculations)

        links_dict[fw_setup_folders] = [fw_chunk_calculations] 
        # FireWork: update database, 
        # (includes reading relaxed structure and energy)
        fw_update_converged_data = update_converged_data(chunk_size = -1)
        workflow_list.append(fw_update_converged_data)
        links_dict[fw_chunk_calculations] =[fw_update_converged_data]

        # eliminate adsorbate pairs too close
        # early exit here
        fw_eliminate_pairs = eliminate_pairs(adsorbate_name = adsorbate_name, bond_length = bond_length)
        workflow_list.append(fw_eliminate_pairs)
        links_dict[fw_update_converged_data] = [fw_eliminate_pairs]

    ### loop ends ###
    wf = Workflow(workflow_list, links_dict)
    return wf


