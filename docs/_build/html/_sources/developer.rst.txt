.. _developer:

For Future Developers
=====================

Become a developer if you want to add more functionality to critcatworks. Raising an issue on github would be the first step. Do not hesitate to contact us if you have any questions.

The goal of critcatworks is to automate nanocluster-surface related research. In critcatworks belongs everything which joins already available building blocks into a complex workflow.

If you want to improve the workflow manager side, Fireworks is the dependency to work on.

If you want to instead create a nanocluster tool, it most likely belongs in cluskit (unless very simple). Since cluskit is also developed in this group, a concerted effort to make that tool available both in cluskit and in critcatworks can be tackled (contact us in that case). Make sure that the functionality is not already in cluskit!


How to Write Custom Firetasks
-----------------------------

Writing your custom Firetask is easy. You just need to wrap your function in a class with some decorations beforehand and afterwards. Before you start implementing your first Firetask, make sure to have a basic knowledge about Fireworks.

.. code-block:: python

    from fireworks import explicit_serialize, FiretaskBase, FWAction

    @explicit_serialize
    class MyCustomTask(FiretaskBase):
        """ 
        Custom Firetask template.

        Args:
            required_parameter1 (any):  you can read any parameters 
                                        during creationg of this task

            required_parameter2 (any):  lists, dictionaries, arrays, etc. are all fine, 
                                        but no pure python objects

            optional_parameter1 (any):  Remember to add them to the list below
        Returns:
            FWAction : Firework action, updates fw_spec
        """
        _fw_name = 'MyCustomTask'
        required_params = ['required_parameter1', 'required_parameter2']
        optional_params = ['optional_parameter1']

        def run_task(self, fw_spec):
            # those values cannot be modified during runtime of the workflow
            optional = self.get("optional_parameter1", "default_value")
            important_parameter = self["required_parameter1"]
            another_parameter = self["required_parameter2"]


            # you can also get information from the firework spec (this can be 
            #modified during runtime of the workflow)
            analysis_ids = fw_spec.get("temp", {}).get("analysis_ids", [1, 2, 3])
            # analysis_ids becomes calc_ids and is stored later
            calc_ids = analysis_ids
            
            # run your custom code
            mycustom_dct = {1 :2, 3 : 4}

            
            # check where this file gets written
            with open('mycustomfile.txt', 'w') as outfile:
                json.dump(mycustom_dct, outfile)

            # fireworks
            # Store information for future jobs to fetch and/or to keep record
            fw_spec["calc_ids"] = calc_ids

            # important to remove those, otherwise they would 
            # overwrite the next Firework's _category and name
            fw_spec.pop("_category")
            fw_spec.pop("name")

            # always return a FWAction object. 
            # other arguments can deviate or defuse the workflow
            return FWAction(update_spec=update_spec)


Reading from and Writing to Permanent database
----------------------------------------------

For interacting with an external database, consider using the functions in :code:`critcatworks.database.extdb`.

The function *get_external_database* connects you to a database using *extdb_connect*.

Then, for instance *fetch_simulations* can get multiple simulations by id.

Lastly, *update_simulations_collection* uploads one simulation document to the database.

For other functionalities in :code:`critcatworks.database.extdb` consult the code documentation.


Fireworks Spec Entries
-----------------------

The current workflows use the following *fw_spec* entries. It is recommended to adhere to the structure but is not prohibited in any way.

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

    calc_paths (list of str)
        paths to the dft calculations, sorted by adsorbate ids
    calc_ids (list of int) 
        ids of simulations in permanent database
    is_converged_list (list of int)
        1 - converged, 
        0 - not converged calculation, 
        same order as calc_paths
    fps_ranking (list of int)
        adsorbate ids ordered by FPS ranking
    analysis_ids (list of int)
        calculation ids which have been analysed and where analysis can be processed

    calc_analysis_ids_dict (dict)
        keys are calculation ids before DFT
        values are calculation ids which have been analysed
    cohesive_energy_dct (dict)
        for each chemical formula key, the value corresponds
        to a dict of simulation indices and cohesive energies
        (total energies if no atomic energies were given)

    descmatrix (str)
        path to numpy array. 2D-matrix descriptor, row representing datapoint
    property (list of str)
        property of interest to machine learning
    last_machine_learning_id (int)
        id of last machine learning step

    reference_energy (float) 
        reference energy for the adsorbate. Can be the total energy of the isolated adsorbate molecule or a different reference point

    free_energy_correction (float) 
        constant shift in free
        energy. This is relevant for the coverage ladder
        target energy range.

    branch_dct (dict)
        keys of parent simulations
        with values being lists of child simulations
    direction (bool)
        1 - adding adsorbate
        0 - removing adsorbate

    ne_dct (dict) 
        stores total energies of all calculations with respect to the number of adsorbates and their ids
    n_adsorbates_root (int) 
        number of adsorbates of the root structure
    n_adsorbates (int) 
        number of adsorbates of the current step
    
    is_return (bool)
        current state of the coverage ladder workflow. If True, 
        the ladder search is on the way back to the root level
    is_new_root (bool)
        If True, the last simulation has resulted in a new
        root simulation

    open_branches (list)
        each element is a tuple containing parent simulation ids and direction 

    root_history (list)
        ordered ids of root simulations during the course of
        the workflow, starting with the start_id

    step_history (list)
        each entry is a tuple of
        a list of calculation ids
        and a direction indicator

    calc_parents (dict)
        keys of simulation ids
        with values being parent simulation ids

    start_id (int)
        unique identifier of the simulation which is used to start the workflow
        
