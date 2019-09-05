Singlesites Workflow
====================

The image below gives an overview how the singlesites
workflow can be applied. This workflow requires
a set of relaxed nanoclusters along with their total energies.
Often times, the nanoclusters workflow is run beforehand. The
singlesites workflow starts after the first stop of the flowchart, the human inpection
and selection. This example of human interference can be automated with ease in the future.


.. image:: ../images/singlesites.svg
   :width: 800


The purpose is to automate mapping adsorption energies on the surface of an arbitrary number of nanoclusters simultaneously. By analysing similar nanoclusters (same elements, similar composition, slightly different shapes) machine-learning can make that effort cheaper.
However, if your nanoclusters are way different, it might be better to split them into different singlesites workflow runs.

What does the singlesites workflow do internally? Here is a short overview.


First, the adsorption sites need to be mapped out. Given an arbitrary nanocluster shape, e.g. Platinum:


.. image:: ../images/pt55.png 
   :width: 600


The surface is triangulated by the aid of the Delaunay algorithm (cluskit functionality). In short, Delaunay divides the volume into non-overlapping tetrahedra, of which the outermost triangles are found. The algorithm is quite robust but can fail in rare cases. There is one free parameter, the maximum length of a triangle side, which can be tweaked to repair the failed cases.


.. image:: ../images/pt.tetrahedra.png
   :width: 600

Besides finding the surface atoms, the triangles define top, bridge and hollow sites. The parameter **adsite_types** takes a list of any combination of "top", "bridge" and/or "hollow" to determine which sites to populate.
A 3-fold-hollow site is the geometrical center of the triangle, a bridge site is the center of a triangle side, and a top site is equal to a triangle's vertex.

.. image:: ../images/pt.porcupine.png 
   :width: 600

The vectors pointing outward are computed by the normals of the triangular faces wheras the directing sign is pre-determined. Hollow site vectors are constructed using its only triangle normal vector. Bridge sites are assigned the average of the vecotors of the two adjacent triangular faces. Top sites finally receive the average vector over all triangles sharing the corresponding top site vertex.

The sites can then be populated, by the adsorbate atom (**adsorbate_name**) at a sensible distance (for adsorbate molecules, take a look at the molsinglesites workflow, see below).

As an unsupervised learning step, all mapped-out sites are compared to each other and their similarity is measured (in feature space of a descriptor (default: SOAP)). The sites are ranked, while the most dissimilar sites come first and (almost) symmetrical equivalents come last. This step will make the machine-learning step much more efficient.

With this order in store, the adsorption energy of those sites is computed in chunks (**chunk_size**).
After every chunk, the workflow checks the machine-learning accuracy (in particular MAE - mean absolute error). If it satisfies a convergence threshold (**threshold**), the workflow stops early, otherwise it continues computing DFT calculations in chunks until (**max_calculations**) have been calculated.

At one glance, these are the singlesite-workflow-specific parameters:

:adsorbate_name (str):  element symbold of the adsorbed atom
:chunk_size (int): 
    number of calculations to be run simulataneously. 
    If -1, all calculations are run at once. Not recommended here.
:max_calculations (int): 
    maximum number of DFT calculations in the workflow
:adsite_types (list):   
    adsorption site types, can contain any combination of
    "top", "bridge", "hollow"
:threshold (float):     
    ML accuracy of convergence criterion. When below, the workflow is defused. 



If you want to place molecular adsorbates on the surface, use a similar workflow called *get_molsinglesites_workflow*

(:meth:`critcatworks.workflows.molsinglesites.get_molsinglesites_workflow`).

The only difference is that you define an atoms object instead of **adsorbate_name**. Direction and distance will be predetermined by the user via the dummy atom "X". Bidentate adsorption is not yet supported.

:adsorbate (dict) :  
    adsorbed molecule as atoms dict. Contains an "X" dummy atom
    which indicates the anchor point to the nanocluster


Define the molecule with ase, it will get converted automatically into a dictionary, e.g.:

.. code-block:: python

    import ase
    pos = np.array([[ 0.00000000e+00,  0.00000000e+00,  1.16489000e-01],
       [ 0.00000000e+00,  9.39731000e-01, -2.71808000e-01],
       [ 8.13831000e-01, -4.69865000e-01, -2.71808000e-01],
       [-8.13831000e-01, -4.69865000e-01, -2.71808000e-01],
       [ 0.00000000e+00, -1.54520895e-06,  1.91648900e+00]])

    adsorbate_x = ase.Atoms('NH3X', positions=pos)
