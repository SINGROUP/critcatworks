from fireworks import LaunchPad, Workflow
import pathlib
import os,time, sys
import logging
import ase
from scipy.spatial.distance import pdist
import getpass
import numpy as np
# internal modules
from critcatworks.workflows import get_molsinglesites_workflow
from critcatworks.database import mylaunchpad

def read_structures_locally(path):
    structures = []
    path = pathlib.Path(path).resolve()
    for idx, p in enumerate(pathlib.Path(path).iterdir()):
        if p.is_file():
            logging.debug("nanocluster path " + str(p) + " stem " + str(p.stem))
            try:
                atoms = ase.io.read(str(p))
                # set cell to 2.5 the diameter
                pos = atoms.get_positions()
                pdist(pos)
                diameter = pdist(pos).max()
                mpl = 2.5

                atoms.set_cell([diameter * mpl, diameter * mpl, diameter * mpl])
                structures.append(atoms)
                logging.debug(atoms)
            except ValueError:
                logging.warning("WARNING: file type not understood" + str(p) )
                continue
            except:
                logging.error("Unexpected error:", sys.exc_info()[0])
    return structures

if __name__ == "__main__":
    IS_QUEUE = True
    USERNAME = "mjcritcat"
    #PASSWORD = getpass.getpass()
    PASSWORD = "heterogeniuscatalysis"

    if IS_QUEUE:
        logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logdir = str(pathlib.Path(".").resolve())
        logging.basicConfig(filename = logdir + "/singlesites_workflow.log", level=logging.INFO)

    # set up the LaunchPad and reset it
    launchpad = mylaunchpad.create_launchpad(USERNAME, PASSWORD, server = "atlas")
    launchpad.reset('', require_password=False)

    structures = read_structures_locally("../nc_structures")
    structures = read_structures_locally("../selected_ptcu_structures")

    # setup nh3 molecule with anchor x
    pos = np.array([[ 0.00000000e+00,  0.00000000e+00,  1.16489000e-01],
       [ 0.00000000e+00,  9.39731000e-01, -2.71808000e-01],
       [ 8.13831000e-01, -4.69865000e-01, -2.71808000e-01],
       [-8.13831000e-01, -4.69865000e-01, -2.71808000e-01],
       [ 0.00000000e+00, -1.54520895e-06,  1.91648900e+00]])

    adsorbate_x = ase.Atoms('NH3X', positions=pos)

    wf = get_molsinglesites_workflow(username = "mjcritcat", 
        password = PASSWORD,
        extdb_ids = [1922,1923,1924,1925,1926,1927],
        template_path = str(pathlib.Path("../templates/cp2k_mm_energy.inp").resolve()), 
        worker_target_path = "/wrk/jagermar/DONOTREMOVE/workflow_runs/molsinglesites/testruns/molsinglesites/selected_ptcu_structures/",
        #structures = structures,
        reference_energy = -1.16195386047558 * 0.5,
        adsorbate = adsorbate_x,
        chunk_size = 10,
        max_calculations = 25,
        adsite_types = ["top", "bridge", "hollow"],
        n_max_restarts = 1,
        skip_dft = False,
        )

    # store workflow 
    launchpad.add_wf(wf)
