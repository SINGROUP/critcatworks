Singlesites Workflow
====================

The image below gives an overview how this
workflow can be applied. This workflow requires
a set of relaxed nanoclusters along with their total energies.
Often times, the nanoclusters workflow is run beforehand. The
singlesites workflow starts after the first stop, the human inpection
and selection. This stop and human interference can be automated 
in the future.


.. image:: ../images/singlesites.svg
   :width: 800


Given an arbitrary nanocluster shape


.. image:: ../images/pt55.png 
   :width: 800


The surface is triangulated by the aid of the delaunay algorithm


.. image:: ../images/pt.tetrahedra.png
   :width: 800


Top, bridge and hollow sites are mapped out.


.. image:: ../images/pt.porcupine.png 
   :width: 800

Lastly, they are compared with each other and their similarity to each
other is measured (in feature space of a descriptor)


