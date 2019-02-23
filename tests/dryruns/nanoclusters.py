from fireworks import LaunchPad, Workflow
import pathlib
import os,time, sys
import logging
import ase
from scipy.spatial.distance import pdist
import getpass

# internal modules
from critcatworks.workflows import get_nanoclusters_workflow
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
    PASSWORD = getpass.getpass()
    if IS_QUEUE:
        logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logdir = str(pathlib.Path(".").resolve())
        logging.basicConfig(filename = logdir + "/nanocluster_workflow.log", level=logging.INFO)

    # set up the LaunchPad and reset it
    launchpad = mylaunchpad.create_launchpad(USERNAME, PASSWORD, server = "serenity")
    launchpad.reset('', require_password=False)

    structures = read_structures_locally("../nc_structures")
    wf = get_nanoclusters_workflow(username = "mjcritcat", password = PASSWORD,
        source_path = None,
        template_path = str(pathlib.Path("../templates/cp2k_mm_energy.inp").resolve()), 
        #worker_target_path = "/wrk/jagermar/DONOTREMOVE/workflow_runs/nanoclusters/testruns/nanoclusters",
        worker_target_path = "../dummy_db/output/",
        structures = structures,
        extdb_ids = None,
        skip_dft = True,
        )

    # store workflow 
    launchpad.add_wf(wf)
