from fireworks import LaunchPad, Workflow
import pathlib
import os,time

# internal modules
from critcatworks.database import start_from_structures, start_from_database, read_structures
from critcatworks.structure import compare_nanoclusters, update_converged_data
from critcatworks.dft import setup_folders, chunk_calculations
from critcatworks.database.format import ase_to_atoms_dict
from critcatworks.database.update import initialize_workflow_data

def get_nanoclusters_workflow(template_path, username, password, 
    worker_target_path = None, structures = None, 
    extdb_ids = None, source_path = None, reference_energy=0.0, 
    atomic_energies = {},
    n_max_restarts = 1, skip_dft = False, extdb_connect = {}):
    """
    Workflow to relax the structure of a set of
    nanoclusters using a simulation software (e.g. CP2K).
    The cohesive energies are calculated and summarized.
    
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
        atomic_energies (dict)  :   used for computing cohesive energies, not required
        n_max_restarts (int)  : number of times the calculation is restarted upon failure
        skip_dft (bool) :   If set to true, the simulation step is skipped in all
                            following simulation runs. Instead the structure is returned unchanged.
        extdb_connect (dict):   dictionary containing the keys host,
                                username, password, authsource and db_name.
        
    Returns:
        fireworks.Workflow : nanocluster Fireworks Workflow object
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
        "simulation_method" : "cp2k",
        "n_max_restarts" : 1,
        "workflow_type" : "relax_nanoclusters",
        "atomic_energies" : atomic_energies,
        }
    fw_init = initialize_workflow_data(username, password, parameters, 
        name = "UNNAMED", workflow_type = "relax_nanoclusters",
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

    # Firework: setup folders for DFT calculations
    fw_setup_folders = setup_folders(target_path = worker_target_path, name = "cp2k_nanoclusters_id")

    # FireWork: setup, run and extract DFT calculation
    # (involves checking for errors in DFT and rerunning)
    fw_chunk_calculations = chunk_calculations(template = template, target_path = worker_target_path, 
        n_max_restarts = n_max_restarts, simulation_method = "cp2k", skip_dft = skip_dft)

    # FireWork: compare stability of nanoclusters. 
    #Computes cohesive energies if atomic_energies of all involved elements are given
    fw_compare_nanoclusters =  compare_nanoclusters(atomic_energies = atomic_energies)

    # add above Fireworks with links
    workflow_list = [fw_init,
        fw_get_structures, 
        fw_setup_folders,
        fw_chunk_calculations,
        fw_compare_nanoclusters,
        ]

    links_dict = {
            fw_init: [fw_get_structures],
            fw_get_structures: [fw_setup_folders],
            fw_setup_folders : [fw_chunk_calculations],
            fw_chunk_calculations : [fw_compare_nanoclusters],
            }

    wf = Workflow(workflow_list, links_dict)
    return wf
