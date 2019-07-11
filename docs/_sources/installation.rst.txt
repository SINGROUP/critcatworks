Installation
============

This will install all dependencies of critcatworks, including fireworks, a cp2k-parser, cluskit, a.o.
For setting up your own mongodb database, please consult 
`www.mongodb.com <https://www.mongodb.com/>`_ .

0. 
virtual environment recommended

.. code-block:: bash

    virtualenv --python=python3 critcatenv`
    source critcatenv/bin/activate


1. 
Download this repository

2. Go inside this repository
:code:`cd critcatworks`

3.  
slightly modified NOMAD CP2K parser Installation

.. code-block:: bash

    pip install numpy cython
    cd nomad_parser/python-common
    pip install -r requirements.txt
    pip install -e .

    cd ../parser-cp2k
    pip install -e .

4. 
critcatworks installation

.. code-block:: bash

    cd ../../
    python3 setup.py install


5. Deactivate virtual environment
:code:`deactivate`


Steps 3 and 4 are automatically executed by 
:code:`./easily_install_all.sh`
