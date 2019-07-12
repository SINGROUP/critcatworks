Simple Coverage Workflow
========================

The simple coverage workflow obtains a stable, simple coverage estimate. It does so by iteratively relaxing the coverage structure by DFT and removing too close adsorbates. So far, only atomic adsorbates are supported.

.. image:: ../images/simplecoverage.svg
  :width: 800


Given a nanocluster structure, the adsorption sites are mapped out as described in the :doc:`singlesites_workflow` tutorial. Similarly **adsorbate_name** defines which atomic adsorbate is placed onto sites classified as "top", "bridge" and/or "hollow" in **adsite_types**.

Since the cluster is would likely be overpopulated, the amount of adsorbates is pre-reduced: either the user defines the target number of remaining adsorbates (**n_remaining**) or a inter-adsorbate minimum distance (**bond_length**). If **n_remaining** is given, **bond_length** is ignored in this step.


After each DFT relaxation, it is checked if any two adsorbate atoms have formed a molecule. Too close adsorbates (specified by **bond_length**) are consequently removed.

If no adsorbates are too close, the search has converged and the workflow stops. Otherwise, it continues for **max_iterations**.


At one glance, these are the simple-coverage-workflow-specific parameters:

:adsorbate_name (str): element symbol of the adsorbed atom
:max_iterations (int): maximum number of iterations in the workflow
:adsite_types (list):   
    adsorption site types, can contain any combination of
    "top", "bridge", "hollow"
:bond_length (float):   
    distance in angstrom under which two adsorbed atoms are 
    considered bound, hence too close
:n_remaining (int): 
    number of adsorbates which should remain after the
    first pre-DFT pruning of the adsorbate coverage

    