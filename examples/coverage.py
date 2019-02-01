from fireworks import LaunchPad, Workflow

import pathlib
import os,time

# internal modules
import mylaunchpad
from critcatworks.clusgeo import get_adsites, eliminate_pairs 
from critcatworks.database import read_structures, update_coverage_data
from critcatworks.dft import setup_coverage_folders, setup_coverage_cp2k

if __name__ == "__main__":
    import logging
    IS_QUEUE = True
    if IS_QUEUE:
        logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logdir = str(pathlib.Path(".").resolve())
        logging.basicConfig(filename = logdir + "/logfile_ranked_adsites.log", level=logging.INFO)

    # set up the LaunchPad and reset it
    launchpad = mylaunchpad.create_launchpad()
    launchpad.reset('', require_password=False)

    wf = reduce_coverage_workflow(
        source_path = str(pathlib.Path("../../tests/dummy_db/nc_structures/").resolve()),
        template_path = str(pathlib.Path("../../tests/dummy_db/templates/coverage_cheap_gopt.inp").resolve()), 
        target_path = str(pathlib.Path("../../tests/dummy_db/output/").resolve()),
        reference_energy = -1.16195386047558 * 0.5,
        adsorbate_name = "H",
        max_iterations = 4,
        bond_length = 1.5,
    )

    # store workflow on launchpad
    launchpad.add_wf(wf)