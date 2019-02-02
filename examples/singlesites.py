from fireworks import LaunchPad, Workflow

import pathlib
import os,time

# internal modules
from critcatworks.workflows import get_singlesites_workflow
from critcat.database import mylaunchpad


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

    wf = get_singlesites_workflow(
        source_path = str(pathlib.Path("../../tests/dummy_db/nc_structures/").resolve()),
        template_path = str(pathlib.Path("../../tests/dummy_db/templates/rankedadsites_cheap_gopt.inp").resolve()), 
        target_path = str(pathlib.Path("../../tests/dummy_db/output/").resolve()),
        reference_energy = -1.16195386047558 * 0.5,
        adsorbate_name = "H",
        chunk_size = 5,
        max_calculations = 15,
        )

    # store workflow 
    launchpad.add_wf(wf)
