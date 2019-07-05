.. _developer:

For Future Developers
==============

Become a developer
Developer tutorial coming soon


Fireworks specs entries
-----------------------

update dict for recent relevant entries
 
```python
fw_spec = {
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
    
```


