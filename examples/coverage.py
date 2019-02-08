from fireworks import LaunchPad, Workflow
import pathlib
import os,time, sys
import logging
import ase
from scipy.spatial.distance import pdist

# internal modules
from critcat.database import mylaunchpad
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
    if IS_QUEUE:
        logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logdir = str(pathlib.Path(".").resolve())
        logging.basicConfig(filename = logdir + "/coverage_workflow.log", level=logging.INFO)

    # set up the LaunchPad and reset it
    launchpad = mylaunchpad.create_launchpad()
    launchpad.reset('', require_password=False)

    wf = get_coverage_workflow(
        source_path = str(pathlib.Path("./nc_structures/").resolve()),
        template_path = str(pathlib.Path("./templates/gopt.inp").resolve()), 
        worker_target_path = "../tests/dummy_db/output/",
        reference_energy = -1.16195386047558 * 0.5,
        adsorbate_name = "H",
        max_iterations = 4,
        adsite_types = ["top"],
        username = "marc",
        n_max_restarts = 1,
        skipt_dft = True,
        bond_length = 1.5,
    )

    # store workflow on launchpad
    launchpad.add_wf(wf)