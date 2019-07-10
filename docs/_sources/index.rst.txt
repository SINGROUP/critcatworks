.. critcatworks documentation master file, created by
   sphinx-quickstart on Mon Jul  1 13:35:01 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Critcatworks
=======================

Workflow manager for DFT simulations on nanocluster databases using Fireworks

This site is UNDER CONSTRUCTION

The project aims to partially automate nanocatalyst discovery. During the course
of this project, critcatworks emerged. It has made screening of catalyst search space
and management and analysis of large datasets easier. The following sketch shows an
overview of the envisioned and partly realized catalyst discovery platform.

.. image:: ./images/critcat_2_flowchart.svg 
   :width: 800


The package features various workflows:

1. singlesites workflow

   a) automatically generate adsorbates on nanoclusters
   b) rank the cluster-adsorbate structures based on similarity
   c) run CP2K calculations
   d) do machine learning on the fly

2. nanocluster workflow

   a) relax nanocluster structures
   b) sort them with respect to stability 

3. coverage workflow

   a) automatically cover nanoclusters with adsorbates
   b) reduce the coverage with a simple heuristic step by step with DFT

4. coverage ladder workflow

   a) starting from a coverage with adsorbates
   b) sophisticated coverage ladder algorithm searches for the optimal adsorbate coverage


Take a look at the :doc:`installation` and :doc:`quickstart` tutorials for the first steps. 


Tutorials
=========

In order to use the workflow package confidently, it is advised to become familiar with
Fireworks: `Documentation <https://materialsproject.github.io/fireworks/>`_.

.. toctree::
   
   installation
   quickstart
   checklist
   tutorials/workflows
   database
   developer   


Contributing
============

Contribute by raising an issue on github if you encounter a problem, or develop your own
nanocluster workflow. Follow the guidelines in the developer tutorial: :ref:`developer`

Code Documentation
==================

.. toctree::

   src/doc/modules


About
=====

Authors
-------
This package is developed at Aalto University, Department of Applied Physics
by the `Surfaces and Interfaces at the Nanoscale (SIN) 
<https://www.aalto.fi/department-of-applied-physics/surfaces-and-interfaces-at-the-nanoscale-sin>`_ 
group.

Contact
-------
If you encounter issues with the software, or want to suggest improvements,
please use the `issues feature <https://github.com/SINGROUP/critcatworks/issues>`_ of
github. 

.. raw:: html

    <script type="text/javascript" src="_static/js/spamspan.js"></script>
    <p>If you want to contact the authors on other matters, use the following email:
        <span class="spamspan">
        <span class="u">marc.jager</span>
        [at]
        <span class="d">aalto [dot] fi</span>
        </span>
    </p>

License
-------
The software is licensed under the GNU General Public License Version 3.0

Funding
-------
This project has received funding from the Finnish Foundation for 
Technology Promotion and the European Union's Horizon 2020 
research under grant agreement number no. 686053 CRITCAT.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
