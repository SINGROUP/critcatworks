# critcatworks

CritCatWorks is a workflow manager for DFT simulations on nanocluster databases 
built on 
[Fireworks](https://materialsproject.github.io/fireworks/).
The documentation along with tutorials is available at
[https://singroup.github.io/critcatworks/](https://singroup.github.io/critcatworks/)


# Quick example

```python
# This example will only run if you have a connection to a mongodb database
# containing a nanocluster simulation (Execute nanocluster workflow to get
# one). username, password, extdb_ids along with template_path and 
# worker_target_path need to be provided

from fireworks import LaunchPad, Workflow
import pathlib
import os, sys
import getpass
from critcatworks.database import mylaunchpad
from critcatworks.workflows.coverage import get_coverage_workflow

USERNAME = "myusername"
PASSWORD = getpass.getpass()

# set up the LaunchPad and reset it
launchpad = mylaunchpad.create_launchpad(USERNAME, PASSWORD)

# simple coverage workflow
wf = get_coverage_workflow(username = USERNAME, 
    password = PASSWORD,
    template_path = str(pathlib.Path("./templates/triton_gopt.inp").resolve()), 
    worker_target_path = "/path/to/my/calculation/directory/",
    extdb_ids = [245,],                          # unique identifier of your nanocluster structure in your mongodb database
    reference_energy = -1.16195386047558 * 0.5,  # half h2 total energy as reference
    adsorbate_name = "H",
    max_iterations = 4,                          # coverage workflow iterations
    adsite_types = ["top", "hollow", ],          # initial population of the adsorption sites on the nanocluster
    n_max_restarts = 1,                          # one DFT restart upon failure
    skip_dft = False,
    bond_length = 1.5,                           # converges when no adsorbates are as close as the specified bond length
    n_remaining = 80,                            # initial reduction of adsorbates to this number 
    extdb_connect = {"db_name": "ncdb"},         # default database is testdb. use this for production runs
)   

# store workflow on launchpad
launchpad.add_wf(wf)
```

# Installation

critcatworks 

0. virtual environment recommended
virtualenv --python=python3 critcatenv
source critcatenv/bin/activate

1. Download this repository

```sh
git clone git@github.com:SINGROUP/dscribe.git
```

2. Go inside this repository
cd critcatworks

3.  slightly modified NOMAD CP2K parser Installation

This step will become obsolete once the cp2k parser
is available through pip

```sh
pip install numpy cython

cd nomad_parser/python-common
pip install -r requirements.txt
pip install -e .

cd ../parser-cp2k
pip install -e .
```

4. critcatworks installation

```sh
cd ../../
python3 setup.py install
```

--

steps 3 and 4 are automatically executed by 
```sh
./easily_install_all.sh
```

