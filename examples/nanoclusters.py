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
    USERNAME = "myusername"
    PASSWORD = getpass.getpass()
    if IS_QUEUE:
        logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logdir = str(pathlib.Path(".").resolve())
        logging.basicConfig(filename = logdir + "/nanocluster_workflow.log", level=logging.INFO)

    # set up the LaunchPad and reset it
    launchpad = mylaunchpad.create_launchpad(USERNAME, PASSWORD, lpadname = "mjfireworkstriton")
    #launchpad.reset('', require_password=False)

    structures = read_structures_locally("./ptx55")
    wf = get_nanoclusters_workflow(username = "myusername", password = PASSWORD,
        source_path = None,
        template_path = str(pathlib.Path("templates/triton_gopt.inp").resolve()), 
        #worker_target_path = "/wrk/jagermar/DONOTREMOVE/workflow_runs/nanoclusters/production/ptcu_selected_clusters",
        worker_target_path = "/scratch/work/jagerm1/workflow_runs/nanoclusters/production/selected_ptni_clusters",
        structures = structures,
        extdb_ids = None,
        skip_dft = False,
        extdb_connect = {"db_name": "ncdb"},
        )

    # store workflow 
    launchpad.add_wf(wf)
