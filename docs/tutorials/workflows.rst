.. _workflows:

Workflows
=========

The automated workflows facilitate generation and analysis of large datasets of
adsorption on nanoclusters.

At first, a set of nanoclusters can be relaxed by DFT. Clusters with the same composition will
be ranked based on their stability. Human selection of clusters for the next step is currently
required.

After settling on a reduced set of nanoclusters, their adsorption sites are detected and
populated. Two options are possible. Either studying a single adsorbate or a coverage of 
adsorbates on the nanoclusters.

For the former, the single adsorbates are computed by DFT in chunks, smartly pre-selected, until
enough points are acquired in order to infer the rest with machine-learning. The remaining points
are predicted with a specified accuracy.

For the latter, a nanocluster is covered by adsorbates using simple heuristics and an iterative
DFT relaxation followed by a removal of a single adsorbate, one by one. This usually leads to
too low coverages but is a good starting point for the next workflow.
Given a starting coverage, the coverage ladder workflow adds and removes adsorbates one by one,
with a parallel DFT computation of addition/removal candidates, and branches out over time to
consistently find a significantly lower-energy coverage.


This figure provides an overview on how the implemented workflows depend on each other.


.. image:: ../images/implemented_workflows_overview.svg
   :width: 800



Most workflows have the following arguments in common: 


:template_path (str):   
    absolute path to input file for calculations. 
    It works as a template which is later modified by the
    simulation-specific Firework.
:username (str): user who executed the workflow
:password (str): password for user to upload to the database
:worker_target_path (str): absolute path on computing resource. Directory needs to exist
:reference_energy (float):  
    reference energy for the adsorbate. Can be the
    total energy of the isolated adsorbate molecule
    or a different reference point
:n_max_restarts (int): number of times the DFT calculation is restarted upon failure
:skip_dft (bool):   
    If set to true, the simulation step is skipped in all
    following simulation runs. Instead the structure is returned unchanged.
:extdb_connect (dict):   
     dictionary containing the keys host,
     username, password, authsource and db_name 
     of the database to connect to. Defaults to
     a test database for critcat.
     If db_name is set to ncdb, this will upload
     the data to the production database.


The following workflows come with a tutorial (there are a few more simple or similar workflows without a tutorial):

.. toctree::
   :maxdepth: 1

   nanocluster_workflow
   singlesites_workflow
   coverage_workflow
   coverageladder_workflow
