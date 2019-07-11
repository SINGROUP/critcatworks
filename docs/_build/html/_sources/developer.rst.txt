.. _developer:

For Future Developers
=====================

Become a developer
Developer tutorial coming soon


Fireworks specs entries
-----------------------


update dict for recent relevant entries


:simulations (dict): 
    simulation collection entries for this workflow.
    Usually, simulations are not stored here, since large amounts
    of documents would slow the workflow manager down

:workflow (dict): 
    relevant information about this workflow,
    entry for workflow collection

:machine_learning (dict): 
    machine_learning instances of this workflow
    entries for machine_learning collection

:n_calcs_started (int): 
    number of calculations which have already been started

:extdb_connect (dict):
    Connection information to permanent mongodb database containing the keys host, username, password, 
    authsource and db_name.

:temp (dict):
    :calc_paths (list of str): 
        paths to the dft calculations, sorted by adsorbate ids
    :calc_ids (list of int): 
        ids of simulations in permanent database
    :is_converged_list (list of int): 
        1 - converged, 
        0 - not converged calculation, 
        same order as calc_paths
    :fps_ranking (list of int): 
        adsorbate ids ordered by FPS ranking
    :analysis_ids (list of int): 
        calculation ids which have been analysed and where analysis can be processed

    :calc_analysis_ids_dict (dict):
        keys are calculation ids before DFT
        values are calculation ids which have been analysed
    :cohesive_energy_dct (dict):
        for each chemical formula key, the value corresponds
        to a dict of simulation indices and cohesive energies
        (total energies if no atomic energies were given)

    :descmatrix (str): 
        path to numpy array. 2D-matrix descriptor, row representing datapoint
    :property (list of str): 
        property of interest to machine learning
    :last_machine_learning_id (int): 
        id of last machine learning step

    :reference_energy (float): reference energy for the adsorbate. Can be the total energy of the isolated adsorbate molecule or a different reference point

    :free_energy_correction (float): constant shift in free
        energy. This is relevant for the coverage ladder
        target energy range.

    :branch_dct (dict):
        keys of parent simulations
        with values being lists of child simulations
    :direction (bool):
        1 - adding adsorbate
        0 - removing adsorbate

    :ne_dct (dict): stores total energies of all calculations with respect to the number of adsorbates and their ids
    :n_adsorbates_root (int): number of adsorbates of the root structure
    :n_adsorbates (int): number of adsorbates of the current step
    
    :is_return (bool): 
        current state of the coverage ladder workflow. If True, 
        the ladder search is on the way back to the root level
    :is_new_root (bool):
        If True, the last simulation has resulted in a new
        root simulation

    :open_branches (list):  
        each element is a tuple containing parent simulation ids and direction 

    :root_history (list):
        ordered ids of root simulations during the course of
        the workflow, starting with the start_id

    :step_history (list):
        each entry is a tuple of
        a list of calculation ids
        and a direction indicator

    :calc_parents (dict):
        keys of simulation ids
        with values being parent simulation ids

    :start_id (int):  
        unique identifier of the simulation which is used to start the workflow
        



TO DELETE

.. code-block:: JSON

    {
        simulations (dict) : "simulation entries for this workflow",
        workflow (dict) : "relevant information about this workflow",
        machine_learning (dict) : "machine_learning instances of this workflow",
        temp (dict) :
            {
            calc_paths (list of str) : "paths to the dft calculations, sorted by adsorbate ids",
            calc_ids (list of int) : "ids of simulations in permanent database",
            is_converged_list (list of int) : "1 - converged, 0 - not converged calculation, same order as calc_paths",
            fps_ranking (list of int) : "adsorbate ids ordered by FPS ranking",
            analysis_ids (list of int) : "calculation ids which have been analysed and where analysis can be processed",
            descmatrix (str) : "path to numpy array. 2D-matrix descriptor, row representing datapoint",
            property (list of str) : "property of interest to machine learning",
            last_machine_learning_id (int) : "id of last machine learning step",
            }
        n_calcs_started (int) : "number of calculations which have already been started",
    }


