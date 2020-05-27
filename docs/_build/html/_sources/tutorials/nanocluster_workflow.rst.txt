Nanocluster Workflow
====================

This workflow simply relaxes the structure of a set of nanoclusters

If atomic energies (**atomic_energies**) are given, the cohesive energies of the relaxed nanoclusters
are computed. Otherwise, the total energy is used to compare 
nanoclusters of the same composition. 

Currently, stable nanoclusters have to be picked manually once the workflow 
has finished, but this can be automated in the future.


.. image:: ../images/nanoclusters.svg
   :width: 800


This is an example how to use the nanocluster workflow:

.. literalinclude:: ../../examples/nanoclusters.py
   :language: python
   :lines: 1-64

As in most workflows, the initial structures can be read in three ways:

:structures:
    list of ase.Atoms objects
:extdb_ids:
    list of unique identifiers of the simulations collection. The
    simulations in the database need to have the correct form
:source_path:
    absolute path on the computing resource to the directory
    where to read the structures from


There are a few workflow-specific arguments:

:atomic_energies (dict):   used for computing cohesive energies, not required

The other arguments are common to all workflows, such as the username and password for the
database or the path to the DFT template. 

