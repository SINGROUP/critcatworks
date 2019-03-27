from fireworks import LaunchPad, Workflow
import pathlib
import os,time, sys
import logging
import ase
from scipy.spatial.distance import pdist
import getpass

# internal modules
from critcatworks.database import mylaunchpad
from critcatworks.workflows.coverage import get_coverage_workflow

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
    import logging
    IS_QUEUE = True
    USERNAME = "mjcritcat"
    PASSWORD = getpass.getpass()
    if IS_QUEUE:
        logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logdir = str(pathlib.Path(".").resolve())
        logging.basicConfig(filename = logdir + "/coverage_workflow.log", level=logging.INFO)

    # set up the LaunchPad and reset it
    launchpad = mylaunchpad.create_launchpad(USERNAME, PASSWORD, lpadname = "mjfireworkstriton")
    #launchpad.reset('', require_password=False)
    
    wf = get_coverage_workflow(username = "mjcritcat", 
        password = PASSWORD,
        template_path = str(pathlib.Path("./templates/triton_gopt.inp").resolve()), 
        worker_target_path = "/scratch/work/jagerm1/workflow_runs/coverage/production/selected_ptni_clusters",
        extdb_ids = [32, 33],
        reference_energy = -1.16195386047558 * 0.5,
        adsorbate_name = "H",
        max_iterations = 4,
        adsite_types = ["bridge", "hollow", ], 
        n_max_restarts = 1,
        skip_dft = False,
        bond_length = 1.5,
        #n_remaining = 100,
        extdb_connect = {"db_name": "ncdb"},
    )

    # store workflow on launchpad
    launchpad.add_wf(wf)
