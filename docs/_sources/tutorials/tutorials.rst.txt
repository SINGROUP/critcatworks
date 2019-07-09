.. _tutorials:

Tutorials
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

The following workflows are available (alongside a few others):

ADD OVERVIEW FIGURE


In order to use the workflow package confidently, it is advised to become familiar with
Fireworks: `Documentation <https://materialsproject.github.io/fireworks/>`_.


.. toctree::
   :maxdepth: 1

   nanocluster_workflow
   singlesites_workflow
   coverage_workflow
   coverageladder_workflow
   developer   
