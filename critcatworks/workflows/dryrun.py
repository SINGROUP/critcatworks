from fireworks import LaunchPad, Workflow

import pathlib
import os,time

# internal modules
from critcatworks.database import read_structures 
from critcatworks.dft.cp2k import setup_dry_folders, setup_dry_cp2k
import mylaunchpad

def run_cp2k_dry(source_path, template_path, target_path = None):
    """
    Workflow to test if cp2k runs without problems
    """
    # FireWork: Read nanocluster structures and initialise a database
    # object containing set information
    abspath = pathlib.Path(source_path).resolve()

    if target_path == None:
        target_path = os.getcwd()
    else:
        target_path = pathlib.Path(target_path).resolve()

    fw_read_structures = read_structures(abspath)

    # Firework: setup folders for DFT calculations
    fw_setup_folders = setup_dry_folders(target_path = target_path)


    # add above Fireworks with links
    workflow_list = [fw_read_structures, 
        fw_setup_folders,
        ]

    links_dict = {
            fw_read_structures : [fw_setup_folders],
            }

    # FireWork: setup, run and extract DFT calculation
    # (involves checking for errors in DFT and rerunning)
    fw_setup_cp2k = setup_dry_cp2k(template_path = template_path)
    workflow_list.append(fw_setup_cp2k)
    links_dict[fw_setup_folders] = [fw_setup_cp2k]

    wf = Workflow(workflow_list, links_dict)
    return wf


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

    wf =  run_cp2k_dry(
        source_path = str(pathlib.Path("../../tests/dummy_db/nc_structures/").resolve()),
        template_path = str(pathlib.Path("../../tests/dummy_db/templates/nothing.template").resolve()), 
        target_path = str(pathlib.Path("../../tests/dummy_db/output/").resolve()),
        )

    # store workflow and launch it locally, single shot
    launchpad.add_wf(wf)