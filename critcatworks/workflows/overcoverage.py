from fireworks import LaunchPad, Workflow

import pathlib
import os,time

# internal modules
import mylaunchpad
from critcatworks.clusgeo import get_adsites, eliminate_pairs 
from critcatworks.database import read_structures, update_coverage_data
from critcatworks.dft import setup_coverage_folders, setup_coverage_cp2k

def reduce_overcoverage_workflow(source_path, template_path, target_path = None, reference_energy=0.0, 
        adsorbate_name='H', max_iterations = 50, bond_length = 1.0):
    """
    Workflow to determine a stable coverage of a nanocluster with single adsorbate atoms. As a first step, 
    adsorbates are put on top, bridge and hollow sites. Once the structure is relaxed by DFT,
    formed adsorbate molecules (pairs of atoms) are replaced by a single adsorbate.
    The procedure is repeated until no adsorbate molecules form.
    """
    # FireWork: Read nanocluster structure and initialise a database
    # object containing set information
    abspath = pathlib.Path(source_path).resolve()

    if target_path == None:
        target_path = os.getcwd()
    else:
        target_path = pathlib.Path(target_path).resolve()

    fw_read_structures = read_structures(abspath)
    # FireWork: Determine adsites and add to database
    # create structure with overcoverage
    fw_get_adsites = get_adsites(
        reference_energy=0.0, 
        adsorbate_name='H', 
        adsite_types = ["top", "bridge", "hollow"],
        )

    # add above Fireworks with links
    workflow_list = [fw_read_structures, 
        fw_get_adsites, 
        ]

    links_dict = {
            fw_read_structures: [fw_get_adsites], 
            }

    ### loop starts ###
    for i in range(max_iterations):
        # Firework: setup folders for DFT calculations,
        fw_setup_coverage_folders = setup_coverage_folders(target_path = target_path, name = "cp2k_coverage_iter_" + str(i))
        workflow_list.append(fw_setup_coverage_folders)
        if i == 0:
            links_dict[fw_get_adsites] = [fw_setup_coverage_folders]
        else:
            links_dict[fw_eliminate_pairs] = [fw_setup_coverage_folders]


        # FireWork: setup, run and extract DFT calculation
        # (involves checking for errors in DFT and rerunning)
        fw_setup_coverage_cp2k = setup_coverage_cp2k(template_path = template_path, target_path = target_path, name = "cp2k_coverage_iter_" + str(i), n_max_restarts = 0)
        workflow_list.append(fw_setup_coverage_cp2k)

        links_dict[fw_setup_coverage_folders] = [fw_setup_coverage_cp2k] 
        # FireWork: update database, 
        # (includes reading relaxed structure and energy)
        fw_update_coverage_data = update_coverage_data()
        workflow_list.append(fw_update_coverage_data)
        links_dict[fw_setup_coverage_cp2k] =[fw_update_coverage_data]

        # eliminate adsorbate pairs too close
        # early exit here
        fw_eliminate_pairs = eliminate_pairs(adsorbate_name = adsorbate_name, bond_length = bond_length)
        workflow_list.append(fw_eliminate_pairs)
        links_dict[fw_update_coverage_data] = [fw_eliminate_pairs]


    ### loop ends ###

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

    wf = reduce_overcoverage_workflow(
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