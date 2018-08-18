# critcatworks
Worklow manager for DFT simulations on nanocluster databases using Fireworks

Work in progress, the first version will feature a workflow to 
  -automatically generate adsorbates on nanoclusters
  -rank the cluster-adsorbate structures based on similarity
  -run CP2K calculations
  -do machine learning on the fly

Dependencies:
fireworks
clusgeo
sklearn



# Fireworks specs for get_adsites_workflow
```python
fw_spec = {
    fps_ranking (list of int) : "adsorbate ids ordered by FPS ranking",
    ads_structures (list of dict) : "adsorbate structures in ase dicts",
    ids_predicted (list of int) : "ids of datapoints which have been predicted by ML",
    descmatrix (list of list of floats) : "2D-matrix descriptor, row representing datapoint",
    adsorbate_name (str) : "adsorbate atom to be placed on surface of nanocluster",
    mae (float) : "mean absolute error of ML test set",
    n_calcs_started (int) : "number of calculations which have already been started",
    relaxed_structure_dict (dict of int : dict) : "relaxed adsorbate structures, with ids as keys and ase dicts as values",
    predicted_energies (list of float) : "predicted energies of datapoints listed in ids_predicted",
    connect_dict (dict of dict of int : list of int) : "connection between cluster ids as keys and a list of adsorbate ids as values (subdivided into adsite_types)",
    reverse_connect_dict (dict of int : int) : "connection between cluster ids as values and adsorbate ids as keys",
    reference_energy (float) : "reference energy to be used when calculating the adsorption energy",
    adsorbate_energies_dict (dict of int : float) : "total energy of the relaxed adsorbate structure",
    adsorbate_energies_list (list of float) : "total energy of the relaxed adsorbate structure sorted by id",
    best_krr_parameters (dict of str : any) : "best parameters for krr on the last ML iteration",
    nc_ids (list of int) : "nanocluster ids",
    nc_names (list of str) : "nanocluster names",
    nc_structures (list of dict) : "nanocluster structures in ase dicts",
    relaxed_structure_list (list of dict) : "relaxed adsorbate structures as ase dicts sorted by ids",
    nc_atomic_numbers (list of int) : "unique sorted list of atomic numbers present in the dataset",
    adsite_types (list of str) : "top, bridge, hollow sites",
    calc_paths (list of str) : "paths to the dft calculations, sorted by adsorbate ids",
    nc_energies (list of float) : "total energies of the nanoclusters sorted by nanocluster id",
    path (str) : "absolute path to nanocluster structure folder",
    reaction_energies_list (list of float) : "reaction energy (using reference_energy) of the relaxed adsorbate structure sorted by adsorbate id",
}
    
```

 
# Installation
# NOMAD CP2K parser Installation
The code is python 2 and python 3 compatible. First download and install
the nomadcore package:

```sh
git clone https://gitlab.mpcdf.mpg.de/nomad-lab/python-common.git
cd python-common
pip install -r requirements.txt
pip install -e .
```

Then download the metainfo definitions to the same folder where the
'python-common' repository was cloned:

```sh
git clone https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info.git
```

Finally download and install the parser:

```sh
git clone https://gitlab.mpcdf.mpg.de/nomad-lab/parser-cp2k.git
cd parser-cp2k
pip install -e .
```


