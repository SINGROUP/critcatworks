from fireworks import LaunchPad, Workflow
import pathlib
import os,time, sys
import logging
import ase

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
    if IS_QUEUE:
        logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logdir = str(pathlib.Path(".").resolve())
        logging.basicConfig(filename = logdir + "/logfile_ranked_adsites.log", level=logging.INFO)

    # set up the LaunchPad and reset it
    launchpad = mylaunchpad.create_launchpad()
    launchpad.reset('', require_password=False)

    structures = read_structures_locally("./nc_structures")
    wf = get_nanoclusters_workflow(
        source_path = None,
        template_path = str(pathlib.Path("templates/rankedadsites_cheap_gopt.inp").resolve()), 
        target_path = str(pathlib.Path("../tests/dummy_db/output/").resolve()),
        structures = structures,
        extdb_ids = None,
        )

    # store workflow 
    launchpad.add_wf(wf)
