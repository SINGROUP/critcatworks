from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time, pathlib, sys
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import ase, ase.io
import cluskit
import dscribe
import numpy as np
import logging
from scipy.spatial.distance import cdist
from copy import deepcopy
from collections import defaultdict

from critcatworks.database import atoms_dict_to_ase, ase_to_atoms_dict
from critcatworks.database import read_descmatrix, write_descmatrix
from critcatworks.database import adsorbate_pos_to_atoms_lst
from critcatworks.database import join_cluster_adsorbate
from critcatworks.database import split_nanocluster_and_adsorbates
from critcatworks.database.extdb import gather_all_atom_types
from critcatworks.database.extdb import update_simulations_collection
from critcatworks.database.extdb import fetch_simulations
from critcatworks.database.extdb import get_external_database, _query_id_counter_and_increment

from ase.visualize import view

def x2_to_x(points, bondlength = 1.5):
    """Helper function to determine the points closest to all other 
    points. It ranks them removing one by one until no points
    are closer than the specified parameter bondlength.

    Args:
        points (2D ndarray) : Points in n-dimensional space
        bondlength (float) :    criterion up to which points should
                                be removed.

    Returns:
        1D ndarray :    ids of the remaining points, ordered by minimum
                        distance
    """
    dmat = cdist(points, points)
    ids = cluskit.cluster._rank_fps(points, K = None, greedy =False)
    dmat = np.triu(dmat[ids, :][:, ids])
    remaining_ids = np.all((dmat > bondlength) | (dmat == 0), axis =0)
    return ids[remaining_ids == 1]

@explicit_serialize
class AdsorbateEliminationTask(FiretaskBase):
    """Firetask to eliminate too close adsorbate atoms on a covered
    nanocluster. The adsorbates
    can be eliminated either based on a minimum distance
    'bond_length' or the closest ones can be eliminated until
    only n_remaining adsorbates are left.

    The workflow is defused if no adsorbate atoms were removed.

    Args:
        adsorbate_name (str) : element symbol of the adsorbed atom
        bond_length (float) :   distance in angstrom under which two adsorbed atoms are 
                                considered bound, hence too close
        n_remaining (int) : number of adsorbates which should remain after the
                            first pre-DFT pruning of the adsorbate coverage
    
    Returns:
        FWAction :  Firework action, updates fw_spec, possible defuses
                    workflow.
    """

    _fw_name = 'AdsorbateEliminationTask'
    required_params = ['adsorbate_name', ]
    optional_params = ['bond_length', 'n_remaining']

    def run_task(self, fw_spec):
        logging.debug(fw_spec)
        adsorbate_name = self['adsorbate_name']
        bond_length = self.get('bond_length', "")
        n_remaining = self.get('n_remaining', "")

        calc_ids = fw_spec["temp"]["calc_ids"]
        new_calc_ids = []
        update_spec = fw_spec
        simulations = fetch_simulations(fw_spec["extdb_connect"], calc_ids)
        # divide adsorbate and nanocluster
        for idx, calc_id in enumerate(calc_ids):
            source_simulation = simulations[str(calc_id)]
            atoms = atoms_dict_to_ase(source_simulation["atoms"])
            dct = copy.deepcopy(source_simulation)
            dct["source_id"] = calc_id
            dct["inp"] = {}
            dct["inp"]["bond_length"] = bond_length
            dct["inp"]["adsorbate_name"] = adsorbate_name
            dct["adsorbates"] = []
            dct["output"] = {}

            adsorbate_ids = []
            reference_ids = []
            site_class_list = []
            site_ids_list = []

            for adsorbate in source_simulation["adsorbates"]:
                adsorbate_ids.extend(adsorbate["atom_ids"])
                reference_ids.append(adsorbate["reference_id"])
                site_class_list.append(adsorbate.get("site_class",""))
                site_ids_list.append(adsorbate.get("site_ids", []))
            adsorbate_atoms = atoms[np.array(adsorbate_ids, dtype = int)]

            cluster_ids = []
            for nanocluster in source_simulation["nanoclusters"]:
                cluster_ids.extend(nanocluster["atom_ids"])
            cluster_atoms = atoms[np.array(cluster_ids, dtype = int)]

            adsorbate_positions = adsorbate_atoms.get_positions()
            
            if n_remaining:
                remaining_ids = cluskit.cluster._rank_fps(adsorbate_positions, K = n_remaining, greedy =False)
            elif bond_length:
                remaining_ids = x2_to_x(adsorbate_positions, bondlength = bond_length)

            else:
                logging.warning("give either argument bond_length or n_remaining")

            print("bondlength", bond_length)
            n_kept = remaining_ids.shape[0]

            new_atoms = cluster_atoms
            for idx in remaining_ids:
                adsorbate = adsorbate_atoms[idx]
                new_atoms, cluster_ids, adsorbate_ids = join_cluster_adsorbate(new_atoms, adsorbate)
                dct["adsorbates"].append(dict({"atom_ids" : adsorbate_ids, "reference_id" : reference_ids[idx], 
                    "site_class" : site_class_list[idx], "site_ids" : site_ids_list[idx]}))
                # info about surface atoms not there
            n_adsorbates = len(dct["adsorbates"]) - len(adsorbate_atoms)            
            dct["operations"] = [dict({"add_adsorbate" : n_adsorbates})]
            if n_adsorbates == 0:
                is_done = True
            else:
                is_done = False
            dct["atoms"] = ase_to_atoms_dict(new_atoms)

            # update simulation dict
            # refresh adsorbates list
            # create new simulation
            simulation = update_simulations_collection(extdb_connect = fw_spec["extdb_connect"], **dct)
            logging.info("simulation after adding single adsorbates")
            logging.info(simulation)

            # update internal workflow data
            simulation_id = simulation["_id"]
            if not is_done:
                new_calc_ids.append(simulation_id)

        update_spec["temp"]["calc_ids"] = new_calc_ids
        update_spec.pop("_category")
        update_spec.pop("name")
        if len(new_calc_ids) == 0:
            defuse_workflow = True
        else:
            defuse_workflow = False

        return FWAction(update_spec=update_spec, defuse_workflow = defuse_workflow)

@explicit_serialize
class PerTypeCoverageCreationTask(FiretaskBase):
    """Task to cover cluster fully with adsorbates based on a type of site.

    Additionally, a reference simulation document is created and
    the descriptors of the sites are stored.

    Args:
        reference_energy (float) : reference energy of the adsorbate
        adsorbate_name  (str) : Adsorbate atom name to be placed
                                on all sites found.
        adsite_types  (list of str) : Can be "top", "bridge" or "hollow".

    Returns:
        FWAction : Firework action, updates fw_spec
    """
    _fw_name = 'PerTypeCoverageCreationTask'
    required_params = ['reference_energy', 'adsorbate_name', 'adsite_types']
    optional_params = []

    def run_task(self, fw_spec):
        adsorbate_name = self["adsorbate_name"]
        adsite_types = self["adsite_types"]
        reference_energy = self["reference_energy"]
        calc_ids = fw_spec["temp"]["calc_ids"]
        simulations = fetch_simulations(fw_spec["extdb_connect"], calc_ids)
        workflow_id = fw_spec.get("workflow", {"_id" : -1 }).get("_id", -1)
        update_spec = fw_spec

        logging.debug(fw_spec)
        desc_lst = []
        new_calc_ids = []

        # create reference of adsorbate in order to store its total energy
        # for later constructing adsorption energies
        reference_simulation = update_simulations_collection(extdb_connect = fw_spec["extdb_connect"], atoms = {}, 
            source_id = -1, workflow_id = workflow_id, 
            nanoclusters = [], adsorbates = [], substrates = [], 
            operations = [""], inp = {"adsorbate_name" : adsorbate_name}, 
            output = {"total_energy" : reference_energy},)
        reference_id = reference_simulation["_id"]

        all_atomtypes = gather_all_atom_types(calc_ids, simulations)

        # looping over nc atoms 
        for idx, calc_id in enumerate(calc_ids):

            ##
            # get source simulation
            source_simulation = copy.deepcopy(simulations[str(calc_id)])
            atoms_dict = source_simulation["atoms"]
            atoms = atoms_dict_to_ase(atoms_dict)
            logging.debug(atoms)

            # running cluskit on cluster
            cluster = cluskit.Cluster(atoms)
            cluster.get_surface_atoms()
            descriptor_setup = dscribe.descriptors.SOAP(species = all_atomtypes, 
                nmax = 9, lmax = 6, rcut=5.0, crossover = True, sparse = False)
            cluster.descriptor_setup = descriptor_setup

            #looping over adsorption site type
            dct = copy.deepcopy(source_simulation)
            dct["output"] = {}
            dct["output"]["surface_atoms"] = {}                
            dct["source_id"] = calc_id
            dct["workflow_id"] = workflow_id
            dct["inp"] = {}
            dct["inp"]["adsite_types"] = adsite_types
            dct["inp"]["adsorbate_name"] = adsorbate_name

            add_adsorbate = 0 
            for adsite_type in adsite_types:
                if adsite_type == "top":
                    adsite_type_int = 1
                elif adsite_type == "bridge":
                    adsite_type_int = 2

                elif adsite_type == "hollow":
                    adsite_type_int = 3
                else:
                    logging.error("adsorption site type unknown, known types are: top, bridge, hollow")
                    exit(1)
                # get adsorption sites for a nanocluster
                adspos = cluster.get_sites(adsite_type_int)
                sites_surface_atoms = cluster.site_surface_atom_ids[adsite_type_int]

                # get descriptor
                desc = cluster.get_sites_descriptor(adsite_type_int)
                for i in range(desc.shape[0]):
                    desc_lst.append(desc[i])

                adsorbate_lst = adsorbate_pos_to_atoms_lst(adspos, adsorbate_name)

                #loop over each adsorbate
                for adsorbate, surface_atoms in zip(adsorbate_lst, sites_surface_atoms):
                    add_adsorbate += 1
                    atoms, _ , adsorbate_ids = join_cluster_adsorbate(atoms, adsorbate)
                    dct["adsorbates"].append(dict({"atom_ids" : adsorbate_ids, "reference_id" : reference_id, 
                        "site_class" : adsite_type_int, "site_ids" : surface_atoms.tolist()}))
                    dct["output"]["surface_atoms"][str(adsorbate_ids)] = surface_atoms.tolist()

            joint_atoms_dict = ase_to_atoms_dict(atoms)
            dct["atoms"] = joint_atoms_dict
            dct["operations"] = [dict({"add_adsorbate" : add_adsorbate})]

            # update external database
            simulation = update_simulations_collection(extdb_connect = fw_spec["extdb_connect"], **dct)
            logging.info("simulation after adsorbate coverage")
            logging.info(simulation)

            # update internal workflow data
            simulation_id = simulation["_id"]
            new_calc_ids.append(simulation_id)
        
        descmatrix = np.array(desc_lst).tolist()
            
        # saves descmatrix as a path to a numpy array
        update_spec["temp"]["descmatrix"] = write_descmatrix(descmatrix)
        update_spec["temp"]["calc_ids"] = new_calc_ids
        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)


@explicit_serialize
class CoverageLadderTask(FiretaskBase):
    """This Task is at the heart of the coverage ladder workflow.
    It manages the following steps:
    - gets lowest energy structures from the previous step
    - checks if a new root structure has been found
    - manages branches
    - computes adsorption free energy
    - decides upon termination based on parameter d
    - decides direction of next step

    Args:
        d (int) : maximum depth of the coverage ladder (termination criterion)
        l (int) : number of low-energy structures to carry over to the next step

    Returns:
        FWAction :  Firework action, updates fw_spec, possible defuses
                    workflow.
    """
    _fw_name = 'CoverageLadderTask'
    required_params = ["l", "d"]
    optional_params = []

    def run_task(self, fw_spec):

        analysis_ids = fw_spec["temp"]["analysis_ids"]
        calc_ids = fw_spec["temp"]["calc_ids"]
        calc_parents = fw_spec["temp"]["calc_parents"]
        branch_dct = fw_spec["temp"]["branch_dct"]
        direction = fw_spec["temp"]["direction"]
        is_return = fw_spec["temp"]["is_return"]
        energies = fw_spec["temp"]["property"]
        ne_dct = fw_spec["temp"]["ne_dct"]
        l = self["l"]
        d = self["d"]
        is_new_root = fw_spec["temp"]["is_new_root"]
        n_adsorbates_root = int(fw_spec["temp"]["n_adsorbates_root"])
        n_adsorbates = int(fw_spec["temp"]["n_adsorbates"])
        open_branches = fw_spec["temp"]["open_branches"]
        root_history = fw_spec["temp"]["root_history"]
        step_history = fw_spec["temp"]["step_history"]
        reference_energy = fw_spec["temp"]["reference_energy "]
        free_energy_correction  = fw_spec["temp"]["free_energy_correction"]

        lowenergy_calc_ids, energies = self.find_lowest_energy_structures(l, calc_ids, energies)

        lowest_energy, idx = np.array(energies).min(), np.argmin(np.array(energies))
        lowest_idx = lowenergy_calc_ids[idx]

        is_new_root = self.check_new_root(is_new_root, branch_dct,
                                       lowest_idx, ne_dct, n_adsorbates_root, n_adsorbates, 
                                       reference_energy = reference_energy, free_energy_correction = free_energy_correction)
        calc_ids = lowenergy_calc_ids

        # action upon new root
        if is_new_root == True:
            calc_ids = [lowest_idx]
            is_return = False
            root_history.append(lowest_idx)
        else:
            # add regular branch to end
            # is_return variable ensures
            # reverse ladder action to same level of root
            if not is_return:
                open_branches.append((tuple(calc_ids), direction))

            # going opposite direction by adding branch to open_branches
            if n_adsorbates != n_adsorbates_root:
                if is_return:
                    open_branches.insert(0, (tuple(calc_ids), direction))
                else:
                    open_branches.insert(0, (tuple(calc_ids), not direction))
                is_return = True
            else:
                is_return = False

        ## Task for direction ##
        if is_new_root == True:
            # decision up or down
            dg_diff, assignment = self.compute_dg_diff(lowest_idx, branch_dct, ne_dct, reference_energy = reference_energy, free_energy_correction = free_energy_correction)
            if dg_diff <= 0:
                direction = 1
            else:
                direction = 0
            n_adsorbates_root = n_adsorbates
            logging.info("NEW ROOT "  + str(calc_ids))

            # add opposite direction to open_branches
            open_branches = [(calc_ids, not direction)]
        else:
            # decision for next branch
            # check for open branches
            direction, calc_ids, n_adsorbates = self.decide_next_branch(open_branches, ne_dct)
            step_history.append((calc_ids, direction))

            # in tuple next_branch
            pass

        fw_spec["temp"]["branch_dct"] = branch_dct
        fw_spec["temp"]["is_return"] = is_return
        fw_spec["temp"]["open_branches"] = open_branches
        fw_spec["temp"]["is_new_root"] = is_new_root
        fw_spec["temp"]["direction"] = direction
        fw_spec["temp"]["calc_ids"] = calc_ids
        fw_spec["temp"]["n_adsorbates_root"] = n_adsorbates_root
        fw_spec["temp"]["n_adsorbates"] = n_adsorbates

        fw_spec.pop("_category")
        fw_spec.pop("name")
        if (abs(n_adsorbates - n_adsorbates_root) == d) and ((n_adsorbates - n_adsorbates_root > 0) == direction):
            defuse_workflow = True
        else:
            defuse_workflow = False
        return FWAction(update_spec=fw_spec, defuse_workflow = defuse_workflow)

    def find_lowest_energy_structures(self, l, calc_ids, energies):
        """Finds the l lowest energy structures based on the 
        given total energies.

        Args:
            l (int) : number of lowest-energy structures
            calc_ids (list) : simulation ids
            energies (list) : total energies of the systems

        Returns:
            tuple : list of l calc_ids, list of l total energies
        """
        ids = np.argsort(energies)[:l] 
        calc_ids = np.array(calc_ids)[ids].tolist()
        energies = np.array(energies)[ids].tolist()
        return calc_ids, energies

    def level_check_new_root(self, idx, ne_dct, n_adsorbates_root, n_adsorbates):
        """Has to compare with root, if it is 
        on the same coverage level.
        The new structure becomes the new root if its total energy
        is lower.
        """
        if n_adsorbates == n_adsorbates_root:
            pass
        else:
            return False
        # requires dictionary of energies
        energies = ne_dct[str(n_adsorbates)]
        is_new_root = energies[str(idx)] <= min(energies.values())
        return is_new_root

    def check_new_root(self, is_new_root, branch_dct, lowest_idx, 
        ne_dct, n_adsorbates_root, n_adsorbates, 
        reference_energy = 0.0, free_energy_correction = 0.0):
        """Checks if current lowest-energy configuration is eligible
        for new root. 
        A special check is done after the first step after root.
        It depends on the direction of the last step and checks the
        sign of the adsorption free energy.
        Otherwise, the method level_check_new_root is called.

        Args:
            is_new_root (bool) :    if True, the current step is the first step after
                                    root
            branch_dct (dict) : parent simulation : list of child simulations
            lowest_idx (int) : lowest-energy simulation id
            ne_dct (dict) : stores total energies of all calculations with respect to
                            the number of adsorbates and their ids
            n_adsorbates_root (int) : number of adsorbates of the root structure
            n_adsorbates (int) : number of adsorbates of the current step
            reference_energy (float) :  reference energy for the adsorbate. Can be the
                                        total energy of the isolated adsorbate molecule
                                        or a different reference point
            free_energy_correction (float) :    free energy correction of the adsorption 
                                                reaction at hand

            Returns:
            bool : new root
        """
        # first iteration after new root has been set
        if is_new_root:
            dg_diff, assignment = self.compute_dg_diff(lowest_idx, branch_dct, ne_dct, reference_energy = reference_energy, free_energy_correction = free_energy_correction )
            if (
                    ((dg_diff <= 0.0) and (assignment == "child")) or
                    ((dg_diff > 0.0) and (assignment == "parent"))
            ):
                # get as new root
                is_new_root = True
            else:
                is_new_root = False
        else:
            # decision for new root - only if number of adsorbates matches with old root!
            is_new_root = self.level_check_new_root(lowest_idx, ne_dct, n_adsorbates_root, n_adsorbates)
        return is_new_root

    def compute_dg_diff(self, idx, branch_dct, ne_dct, 
        reference_energy = 0.0, free_energy_correction = 0.0):
        """Helper function to compute the adsorption free energy

        In the case of hydrogen the property is:
        dGdiff(nH) = G(nH) - G(nH-1)) - 0.5 G(H2)

        Args:
            idx (int) : simulation id
            branch_dct (dict) : parent simulation : list of child simulations
            ne_dct (dict) : stores total energies of all calculations with respect to
                            the number of adsorbates and their ids
            reference_energy (float) :  reference energy for the adsorbate. Can be the
                                        total energy of the isolated adsorbate molecule
                                        or a different reference point
            free_energy_correction (float) :    free energy correction of the adsorption 
                                                reaction at hand

            Returns:
            tuple : adsorption free energy (float), 
                    assignment of property to 'parent' or 'child' simulation
        """
        # search where idx appears in branch_dct
        for key, ids in zip(branch_dct.keys(), branch_dct.values()):
            if int(idx) in list(ids):
                parent_idx = key
                break
        # get energies     
        for n_adsorbates, v in ne_dct.items():
            if str(idx) in v.keys():
                energy = v[str(idx)]
                child_n_adsorbates = int(n_adsorbates)
            if str(parent_idx) in v.keys():
                parent_energy = v[str(parent_idx)]
                parent_n_adsorbates = int(n_adsorbates)

        assert (parent_n_adsorbates - child_n_adsorbates) == 1 or (parent_n_adsorbates - child_n_adsorbates) == -1

        if parent_n_adsorbates > child_n_adsorbates:
            # dg diff
            dg_diff = parent_energy - energy - reference_energy + free_energy_correction 
            assignment = "parent"
        else:
            dg_diff = energy - parent_energy - reference_energy + free_energy_correction 
            assignment = "child"
        print("#################################################")
        print("dG_diff", dg_diff, assignment)
        print("dG_diff components", energy, parent_energy, reference_energy ,free_energy_correction )
        print("n_adsorbates", "parent", parent_n_adsorbates, "child", child_n_adsorbates)
        print("#################################################")
        return dg_diff, assignment

    def decide_next_branch(self, open_branches, ne_dct):
        """The next branch from the collected possible branches is
        picked. The first element is removed from the list
        which contains information about direction and
        parent simulation ids

        Args:
            open_branches (list) :  each element is a tuple containing 
                                    parent simulation ids and direction
            ne_dct (dict) : stores total energies of all calculations with respect to
                            the number of adsorbates and their ids

        Returns:
            tuple : direction (bool),
                    simulation ids (list),
                    updated number of adsorbates (int)
        """
        next_branch = open_branches.pop(0)
        direction = next_branch[1]
        calc_ids = next_branch[0]
        for n_adsorbates, v in ne_dct.items():
            if str(calc_ids[0]) in v.keys():
                new_n_adsorbates = int(n_adsorbates)
        return direction, calc_ids, new_n_adsorbates


@explicit_serialize
class GatherLadderTask(FiretaskBase):
    """Firetask for Gathering properties for coverage ladder.
    Update database from previous calculations. 
    Computes properties of systems. Currently, only
    reaction energies (adsorption energies) are computed.
    The ids of the input structures in calc_ids are 
    replaced by the ids of the post-simulation structures 
    from analysis_ids.

    Args:

    Returns:
        FWAction : Firework action, updates fw_spec
    """

    _fw_name = 'GatherLadderTask'
    required_params = []
    optional_params = []

    def run_task(self, fw_spec):
        # NOTE also here analysis_ids needs to be in the right order
        calc_analysis_ids_dict = fw_spec["temp"]["calc_analysis_ids_dict"]
        analysis_ids = fw_spec["temp"]["analysis_ids"]
        calc_ids = fw_spec["temp"]["calc_ids"]
        calc_parents = fw_spec["temp"]["calc_parents"]
        branch_dct = fw_spec["temp"]["branch_dct"]
        direction = fw_spec["temp"]["direction"]
        is_return = fw_spec["temp"]["is_return"]
        ne_dct = fw_spec["temp"]["ne_dct"]
        n_adsorbates = fw_spec["temp"]["n_adsorbates"]

        # reorder analysis_ids
        reordered_analysis_ids = []
        for calc_id in calc_ids:
            analysis_id = calc_analysis_ids_dict[str(calc_id)]
            reordered_analysis_ids.append(analysis_id)

        analysis_ids =  reordered_analysis_ids

        simulations = fetch_simulations(fw_spec["extdb_connect"], analysis_ids)
        energies = []
        for idx in analysis_ids:
            simulation = simulations[str(idx)]
            energies.append(simulation["output"]["total_energy"]) 

        # find calc_id in calc_parents and replace with analysis_id
        for _, children in calc_parents.items():
            for i, child in enumerate(children):
                if child in calc_ids:
                    idx = calc_ids.index(child)
                    children[i] = analysis_ids[idx]
        calc_ids = analysis_ids

        # add calc_parents to branch_dct
        if is_return:
            current_direction = not direction
        else:
            current_direction = direction

        for parent, children in calc_parents.items():
            branch_dct.setdefault(parent, []).extend(children)

        # add total energies to dictionary
        for calc_id, energy in zip(calc_ids, energies):
            ne_dct.setdefault(str(n_adsorbates), {})[str(calc_id)] = energy


        fw_spec["temp"]["branch_dct"] = branch_dct
        fw_spec["temp"]["property"] = energies 
        fw_spec["temp"]["calc_ids"] = calc_ids
        fw_spec["temp"]["ne_dct"] = ne_dct
        fw_spec["temp"]["analysis_ids"] = [] 
        fw_spec["temp"]["calc_parents"] = calc_parents
        fw_spec.pop("_category")
        fw_spec.pop("name")
        return FWAction(update_spec=fw_spec)


@explicit_serialize
class AddRemoveLadderTask(FiretaskBase):
    """Firetask in the coverage ladder workflow to either add
    or remove an adsorbate (at k different positions)

    The decision where to add or remove adsorbates is
    quantified by the ranking_metric, either through
    a similarity or spatial distance metric.

    The bond_length parameters ensures that the positions
    are not too cramped, hence adsorbates never 
    come too close to each other.

    Args:
        bond_length (float) :   distance in angstrom under which two adsorbed atoms are 
                                considered bound, hence too close
        k (int) :   number of empty candidate sites for adding / 
                    adsorbed atoms for removing to consider per step
        ranking_metric (str) : 'similarity' or 'proximity'. Metric based on which to choose
                                k candidates (empty sites / adsorbates)

    Returns:
        FWAction : Firework action, updates fw_spec
    """
    _fw_name = 'AddRemoveLadderTask'
    required_params = ["bond_length", "k", "ranking_metric"]
    optional_params = []

    def run_task(self, fw_spec):
        bond_length = self["bond_length"]
        k = self["k"]
        ranking_metric = self["ranking_metric"]
        direction = fw_spec["temp"]["direction"]
        calc_ids = fw_spec["temp"]["calc_ids"]
        n_adsorbates = fw_spec["temp"]["n_adsorbates"]
        workflow_id = fw_spec.get("workflow", {"_id" : -1 }).get("_id", -1)
        db = get_external_database(fw_spec["extdb_connect"])
        calc_parents = {}
        new_calc_ids = []
        simulations = fetch_simulations(fw_spec["extdb_connect"], calc_ids)
        for parent_id, source_simulation  in simulations.items():
            simulation = deepcopy(source_simulation)
            if direction == 1:
                new_ids = self.add_adsorbate(simulation, bond_length, db, ranking_metric = ranking_metric, k = k, workflow_id = workflow_id)
            else:
                new_ids = self.remove_adsorbate(simulation, bond_length, db, ranking_metric = ranking_metric, k = k, workflow_id = workflow_id)
            calc_parents[parent_id] = new_ids
            new_calc_ids.extend(new_ids)

        if direction == 1:
            n_adsorbates += 1
        else:
            n_adsorbates -= 1

        logging.info("number of adsorbates: " + str(n_adsorbates))
        logging.info("after adding / removing one adsorbate: " + str(direction))

        fw_spec["temp"]["calc_parents"] = calc_parents
        fw_spec["temp"]["calc_ids"] = new_calc_ids
        fw_spec["temp"]["n_adsorbates"] = n_adsorbates
        fw_spec.pop("_category")
        fw_spec.pop("name")
        return FWAction(update_spec=fw_spec)

    def add_adsorbate(self, simulation, bond_length, db, 
        ranking_metric = "proximity", k = 7, 
        adsorbate_name = "H", reference_id = -1, workflow_id = -1):
        """Adds one adsorbate to the structure at k different positions,
        predefined by adsorption sites.

        The new structures are added to the simulation collection.

        Only supports atomic adsorbates so far.

        Args:
            simulation (dict) : simulation document in custom format
            bond_length (float) :   distance in angstrom under which two adsorbed atoms are 
                                    considered bound, hence too close
            db (pymongo object) : connection to the mongodb database
            k (int) :   number of empty candidate sites for adding / 
                        adsorbed atoms for removing to consider per step
            ranking_metric (str) : 'similarity' or 'proximity'. Metric based on which to choose
                                    k candidates (empty sites / adsorbates)

            adsorbate_name (str) :  element symbol of the adsorbed atom
            reference_id (int) :    simulation id of the adsorbate reference
            workflow_id (int) :  unique identifier of the current workflow

        Returns:
            list : new simulation ids
        """
        print("adding 1 adsorbate")
        atoms_dict = simulation["atoms"]
        atoms = atoms_dict_to_ase(atoms_dict)
        cluster_atoms, adsorbate_atoms, site_ids_list, site_class_list, reference_ids, adsorbate_ids = split_nanocluster_and_adsorbates(simulation)

        # get sites of nanocluster
        cluster = cluskit.Cluster(cluster_atoms)
        adsorption_sites = cluster.get_sites(-1)

        # find "empty" sites
        ids = self.get_empty_sites(adsorption_sites, adsorbate_atoms.get_positions(), bond_length= bond_length)

        # rank "empty" sites based on
        # a) soap similarity
        if ranking_metric == "similarity":
            desc = cluster.get_sites_descriptor(-1)
            remaining_desc = desc[ids]
            desc_ids = cluskit.cluster._rank_fps(remaining_desc, K=None, greedy=False)
            ids = ids[desc_ids]

        # b) proximity
        else:
            # already ordered
            ids = ids

        # keep first k sites
        if len(ids) < k:
            remaining_ids = ids
        else:
            remaining_ids = ids[:k]

        # construct structures
        adsorbate_lst = adsorbate_pos_to_atoms_lst(adsorption_sites[remaining_ids], adsorbate_name)

        # create new simulations
        # loop over each adsorbate
        simulations_chunk_list = []
        new_calc_ids = []
        for adsorbate in adsorbate_lst:
            # adsites_dict
            joint_atoms, cluster_ids, adsorbate_ids = join_cluster_adsorbate(atoms, adsorbate)
            joint_atoms_dict = ase_to_atoms_dict(joint_atoms)

            # update external database
            dct = deepcopy(simulation)
            # calculation originated from this:
            dct["source_id"] = simulation["_id"]
            dct["workflow_id"] = workflow_id
            dct["atoms"] = joint_atoms_dict
            dct["operations"] = [dict({"add_adsorbate": 1})]
            dct["adsorbates"].append(dict({"atom_ids": adsorbate_ids, "reference_id": reference_id}))
            # empty previous input
            dct["inp"] = {}
            # TODO add info about site 
            #dct["inp"]["adsite_type"] = adsite_type
            dct["inp"]["adsorbate"] = adsorbate_name
            # empty previous output
            dct["output"] = {}
            #dct["output"]["surface_atoms"] = surface_atoms.tolist()

            # getting only id for uploading simulations in chunks
            dct["_id"] = _query_id_counter_and_increment('simulations', db)

            simulations_chunk_list.append(dct)

            # update internal workflow data
            simulation_id = dct["_id"]
            new_calc_ids.append(simulation_id)

        db["simulations"].insert_many(simulations_chunk_list)
        return new_calc_ids

    def remove_adsorbate(self, simulation, bond_length, db, 
        ranking_metric = "similarity", k = 7, adsorbate_name = "H", 
        workflow_id = -1):
        """Remove one adsorbate from the structure at k different positions.

        The new structures are added to the simulation collection.

        Only supports atomic adsorbates so far.

        Args:
            simulation (dict) : simulation document in custom format
            bond_length (float) :   distance in angstrom under which two adsorbed atoms are 
                                    considered bound, hence too close
            db (pymongo object) : connection to the mongodb database
            k (int) :   number of empty candidate sites for adding / 
                        adsorbed atoms for removing to consider per step
            ranking_metric (str) : 'similarity' or 'proximity'. Metric based on which to choose
                                    k candidates (empty sites / adsorbates)

            adsorbate_name (str) :  element symbol of the adsorbed atom
            workflow_id (int) :  unique identifier of the current workflow

        Returns:
            list : new simulation ids
        """
        atoms_dict = simulation["atoms"]
        atoms = atoms_dict_to_ase(atoms_dict)
        cluster_atoms, adsorbate_atoms, site_ids_list, site_class_list, reference_ids, adsorbate_ids = split_nanocluster_and_adsorbates(
            simulation)
        print("removing 1 adsorbate")

        # get descriptor of adsorbate atoms
        cluster = cluskit.Cluster(cluster_atoms + adsorbate_atoms)
        #cluster.get_sites(-1)
        desc = cluster.get_cluster_descriptor()[len(cluster_atoms) :]
        desc.shape

        # rank adsorbates based on
        # a) soap similarity
        # b) proximity
        if ranking_metric == "similarity":
            ids = cluskit.cluster._rank_fps(desc, K=None, greedy=False)
        else:
            ids = x2_to_x(adsorbate_atoms.get_positions(), bond_length = bond_length)

        # keep first k adsorbates
        if len(ids) < k:
            remaining_ids = ids
        else:
            remaining_ids = ids[:k]

        # construct structures
        adsorbate_lst = adsorbate_pos_to_atoms_lst(adsorbate_atoms.get_positions()[remaining_ids], adsorbate_name)

        # create new simulations
        # loop over each adsorbate
        simulations_chunk_list = []
        new_calc_ids = []

        for idx in remaining_ids:
            print("idx of atom to remove", idx)
            joint_atoms = cluster_atoms
            dct = deepcopy(simulation)
            dct["adsorbates"] = []
            for i, adsorbate_id in zip(np.arange(len(adsorbate_ids)), adsorbate_ids):
                if idx == i:
                    # removed atom is excluded
                    pass
                else:
                    adsorbate = adsorbate_atoms[i]
                    joint_atoms, cluster_ids, new_adsorbate_ids = join_cluster_adsorbate(joint_atoms, adsorbate)
                    dct["adsorbates"].append(dict({"atom_ids": new_adsorbate_ids, "reference_id": reference_ids[idx],
                                                   "site_class": site_class_list[idx], "site_ids": site_ids_list[idx]}))
            # info about surface atoms not there

            # adsites_dict
            print("length of structure", joint_atoms)
            joint_atoms_dict = ase_to_atoms_dict(joint_atoms)

            # calculation originated from this:
            dct["source_id"] = simulation["_id"]
            dct["workflow_id"] = workflow_id
            dct["atoms"] = joint_atoms_dict
            dct["operations"] = [dict({"add_adsorbate": -1})]
            # empty previous input
            dct["inp"] = {}
            # TODO add info about site
            #dct["inp"]["adsite_type"] = adsite_type
            #dct["inp"]["adsorbate"] = adsorbate_name
            # empty previous output
            dct["output"] = {}
            #dct["output"]["surface_atoms"] = surface_atoms.tolist()

            # getting only id for uploading simulations in chunks
            dct["_id"] = _query_id_counter_and_increment('simulations', db)

            # simulation = update_simulations_collection(extdb_connect = fw_spec["extdb_connect"], **dct)
            simulations_chunk_list.append(dct)

            # update internal workflow data
            simulation_id = dct["_id"]
            new_calc_ids.append(simulation_id)

        db["simulations"].insert_many(simulations_chunk_list)
        return new_calc_ids

    def get_empty_sites(self, site_positions, 
        adsorbate_positions, bond_length = 1.5):
        """From a set of predefined adsorption sites, it is reduced
        to contain only sites which are not too close to existing
        adsorbate atoms (with respect to parameter bond_length).
        
        Args:
            site_positions (2D ndarray) : adsorption site positions
            adsorbate_positions (2D ndarray) : positions of all adsorbate atoms
            bond_length (float) :   distance criterion up to which a site
                                    is considered occupied by neighbouring atoms
        
        Returns:
            1D ndarray :    ids of the remaining site positions, ordered by 
                            minimum distance to the adsorbates
        """
        dmat = cdist(site_positions, adsorbate_positions)
        mindist = dmat.min(axis = 1)
        ids = np.argsort(-mindist)
        mindist = mindist[ids]
        remaining_ids = mindist > bond_length
        return ids[remaining_ids == 1]

@explicit_serialize
class NewLadderRootTask(FiretaskBase):
    """Task to start CoverageLadder Workflow. It initiates a few
    fw_spec["temp"] entries:
    -"step_history" 
    -"root_history"
    -"is_return"
    -"start_id"
    -"is_new_root"
    -"n_adsorbates"
    -"n_adsorbates_root"
    -"branch_dct"
    -"open_branches" 
    -"ne_dct" 
    -"direction"
    -"calc_ids"
    -"reference_energy"
    -"free_energy_correction"

    Args:
        start_ids (list) :  unique identifiers of the simulations collection which
                            are used to start the workflow
        reference_energy (float) :  reference energy for the adsorbate. Can be the
                                    total energy of the isolated adsorbate molecule
                                    or a different reference point
        free_energy_correction (float) :    free energy correction of the adsorption 
                                            reaction at hand
        initial_direction (bool) :  True will force the initial step to add an adsorbate,
                                    False will force the initial step to remove an adsorbate

    Returns:
        FWAction : Firework action, updates fw_spec    
    """
    _fw_name = 'NewLadderRootTask'
    required_params = ['start_ids', 'reference_energy', 'free_energy_correction']
    optional_params = ['initial_direction']

    def run_task(self, fw_spec):
        start_ids = self["start_ids"]
        initial_direction = self["initial_direction"]
        reference_energy = self["reference_energy"]
        free_energy_correction = self["free_energy_correction"]
        # currently, only one structure per workflow supported
        start_id = start_ids[0]


        simulations = fetch_simulations(fw_spec["extdb_connect"], start_ids)
        start_simulation = simulations[str(start_id)]
        atoms = start_simulation["atoms"]
        n_adsorbates = len(start_simulation["adsorbates"])
        energy = start_simulation["output"].get("total_energy", 0.0)
        if energy == 0.0:
            logging.warning("total_energy for starting structure not found! Will continue without it, however, first two steps might go in the wrong direction")
        n_adsorbates_root = n_adsorbates
        ne_dct = defaultdict(dict) # dictionary of energies, per number of adsorbates
        ne_dct[str(n_adsorbates)][str(start_id)] = energy

        branch_dct = {"root" : np.array(start_ids)}

        # direction to start from root
        direction = initial_direction
        open_branches = [(start_ids, not direction)]
        root_history = [start_id]


        fw_spec["temp"]["step_history"] = [(start_ids, direction)] 
        fw_spec["temp"]["root_history"] = root_history
        fw_spec["temp"]["is_return"] = False,
        fw_spec["temp"]["start_id"] = start_id,
        fw_spec["temp"]["is_new_root"] = True,
        fw_spec["temp"]["n_adsorbates"] = n_adsorbates
        fw_spec["temp"]["n_adsorbates_root"] = n_adsorbates_root
        fw_spec["temp"]["branch_dct"] = branch_dct
        fw_spec["temp"]["open_branches"] = open_branches
        fw_spec["temp"]["ne_dct"] = ne_dct
        fw_spec["temp"]["direction"] = direction
        fw_spec["temp"]["calc_ids"] = start_ids
        fw_spec["temp"]["reference_energy "] = reference_energy 
        fw_spec["temp"]["free_energy_correction"] = free_energy_correction
        fw_spec.pop("_category")
        fw_spec.pop("name")
        return FWAction(update_spec=fw_spec)


def start_coverage_ladder(start_ids, initial_direction = 1, reference_energy = 0.0, free_energy_correction = 0.0):
    """Firework to start CoverageLadder Workflow. It initiates a few
    fw_spec["temp"] entries:
    -"step_history" 
    -"root_history"
    -"is_return"
    -"start_id"
    -"is_new_root"
    -"n_adsorbates"
    -"n_adsorbates_root"
    -"branch_dct"
    -"open_branches" 
    -"ne_dct" 
    -"direction"
    -"calc_ids"
    -"reference_energy"
    -"free_energy_correction"

    Args:
        start_ids (list) :  unique identifiers of the simulations collection which
                            are used to start the workflow
        reference_energy (float) :  reference energy for the adsorbate. Can be the
                                    total energy of the isolated adsorbate molecule
                                    or a different reference point
        free_energy_correction (float) :    free energy correction of the adsorption 
                                            reaction at hand
        initial_direction (bool) :  True will force the initial step to add an adsorbate,
                                    False will force the initial step to remove an adsorbate

    Returns:
        Firework : Firework NewLadderRootWork
    """
    firetask1  = NewLadderRootTask(start_ids = start_ids, initial_direction = initial_direction, 
        reference_energy = reference_energy, free_energy_correction = free_energy_correction)
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'NewLadderRootTask'},
             name = 'NewLadderRootWork')
    return fw

def add_remove_adsorbate(bond_length, k = 7, ranking_metric = "similarity"):
    """Firework in the coverage ladder workflow to either add
    or remove an adsorbate (at k different positions)

    The decision where to add or remove adsorbates is
    quantified by the ranking_metric, either through
    a similarity or spatial distance metric.

    The bond_length parameters ensures that the positions
    are not too cramped, hence adsorbates never 
    come too close to each other.

    Args:
        bond_length (float) :   distance in angstrom under which two adsorbed atoms are 
                                considered bound, hence too close
        k (int) :   number of empty candidate sites for adding / 
                    adsorbed atoms for removing to consider per step
        ranking_metric (str) : 'similarity' or 'proximity'. Metric based on which to choose
                                k candidates (empty sites / adsorbates)

    Returns:
        Firework : Firework AddRemoveLadderWork
    """
    firetask1  = AddRemoveLadderTask(bond_length = bond_length, k = k, ranking_metric = ranking_metric)
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'AddRemoveLadderTask'},
             name = 'AddRemoveLadderWork')    
    return fw

def gather_ladder():
    """Firework for Gathering properties for coverage ladder.
    Update database from previous calculations. 
    Computes properties of systems. Currently, only
    reaction energies (adsorption energies) are computed.
    The ids of the input structures in calc_ids are 
    replaced by the ids of the post-simulation structures 
    from analysis_ids.

    Args:

    Returns:
        Firework : Firework GatherLadderWork
    """
    firetask1  = GatherLadderTask()
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'GatherLadderTask'},
             name = 'GatherLadderWork')    
    return fw

def step_coverage_ladder(l = 2, d = 4):
    """This Firework is at the heart of the coverage ladder workflow.
    It manages the following steps:
    - gets lowest energy structures from the previous step
    - checks if a new root structure has been found
    - manages branches
    - computes adsorption free energy
    - decides upon termination based on parameter d
    - decides direction of next step

    Args:
        d (int) : maximum depth of the coverage ladder (termination criterion)
        l (int) : number of low-energy structures to carry over to the next step

    Returns:
        Firework : Firework CoverageLadderWork
    """
    firetask1  = CoverageLadderTask(l = l, d = d)
    fw = Firework([firetask1], spec = {'_category' : "lightweight", 'name' : 'CoverageLadderTask'},
             name = 'CoverageLadderWork')    
    return fw


def get_per_type_coverage(reference_energy = 0.0, adsorbate_name='H', adsite_types = ["top", "bridge", "hollow"], 
    descriptor = "soap", descriptor_params = {"nmax" : 9, "lmax" :6, "rcut" : 5.0, 
            "crossover" : True, "sparse" : False}):
    """Firework to cover cluster fully with adsorbates based on a type of site.

    Additionally, a reference simulation document is created and
    the descriptors of the sites are stored.

    Args:
        reference_energy (float) : reference energy of the adsorbate
        adsorbate_name  (str) : Adsorbate atom name to be placed
                                on all sites found.
        adsite_types  (list of str) : Can be "top", "bridge" or "hollow".
        descriptor (str) :  type of descriptor to be used. For a list of
                            descriptors, see the documentation of dscribe
                            Defaults to 'soap'
        descriptor_params (dict) : descriptor parameters

    Returns:
        Firework : Firework PerTypeCoverageCreationWork
    """
    firetask1  = PerTypeCoverageCreationTask(
        adsorbate_name=adsorbate_name, 
        adsite_types = adsite_types,
        reference_energy = reference_energy,
        descriptor = descriptor,
        descriptor_params = descriptor_params
        )
    fw = Firework([firetask1], 
        spec={'_category' : "lightweight", 'name' : 'PerTypeCoverageCreationTask'},
        name = 'PerTypeCoverageCreationWork'
        )
    return fw

def eliminate_pairs(adsorbate_name, bond_length):
    """Firework to eliminate too close adsorbate atoms on a covered
    nanocluster. The adsorbates closest to each other are 
    eliminated based on a minimum distance 'bond_length'.

    The workflow is defused if no adsorbate atoms were removed.

    Args:
        adsorbate_name (str) : element symbol of the adsorbed atom
        bond_length (float) :   distance in angstrom under which two adsorbed atoms are 
                                considered bound, hence too close
    
    Returns:
        Firework : Firework AdsorbateEliminationWork
    """
    firetask1  = AdsorbateEliminationTask(adsorbate_name = adsorbate_name, bond_length = bond_length)
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'AdsorbateEliminationTask'},
             name = 'AdsorbateEliminationWork')
    return fw

def eliminate_closest(adsorbate_name, n_remaining):
    """Firework to eliminate too close adsorbate atoms on a covered
    nanocluster. The adsorbates closest to each other are 
    eliminated until only n_remaining adsorbates are left.

    The workflow is defused if no adsorbate atoms were removed.

    Args:
        adsorbate_name (str) : element symbol of the adsorbed atom
        n_remaining (int) : number of adsorbates which should remain after the
                            first pre-DFT pruning of the adsorbate coverage
    
    Returns:
        Firework : Firework AdsorbateEliminationWork
    """
    firetask1  = AdsorbateEliminationTask(adsorbate_name = adsorbate_name, n_remaining = n_remaining)
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'AdsorbateEliminationTask'},
             name = 'AdsorbateEliminationWork')
    return fw
