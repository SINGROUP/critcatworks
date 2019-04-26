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
from critcatworks.database.extdb import update_simulations_collection
from critcatworks.database.extdb import fetch_simulations
from critcatworks.database.extdb import get_external_database, _query_id_counter_and_increment

from ase.visualize import view

def join_cluster_adsorbate(cluster, adsorbate):
    joint_atoms = cluster + adsorbate
    cluster_ids = list(range(len(cluster)))
    adsorbate_ids = list(range(len(cluster_ids), len(joint_atoms)))

    return joint_atoms, cluster_ids, adsorbate_ids

def adsorbate_pos_to_atoms_lst(adspos, adsorbate_name):
    """
    works with only one adsorbate atom.
    TODO generalize cluskit to return list of adsorbates
    already in ase format.
    Makes this function superfluous.
    """
    atoms_lst = []
    ads_structures_dict = []
    for adsorbate in adspos:
        logging.debug(adsorbate_name)
        logging.debug(adsorbate)
        logging.debug(adsorbate.shape)
        atoms = ase.Atoms(symbols=adsorbate_name, positions=adsorbate.reshape((1,3)))
        atoms_lst.append(atoms)
    return atoms_lst

def gather_all_atom_types(calc_ids, simulations):
    # going through nc atoms once to find atom types
    atomic_numbers = []
    for idx, calc_id in enumerate(calc_ids):
        atoms_dict = simulations[str(calc_id)]["atoms"]
        atoms = atoms_dict_to_ase(atoms_dict)
        atomic_numbers.extend(atoms.get_atomic_numbers())

    sorted_list_atomic_numbers = list(sorted(set(atomic_numbers)))

    all_atomtypes = sorted_list_atomic_numbers
    return all_atomtypes

def sparse_subset(points, r):
    """Return a maximal list of elements of points such that no pairs of
    points in the result have distance less than r.
    """
    result = []
    ids = []
    for idx, p in enumerate(points):
        if all(dist(p, q) >= r for q in result):
            result.append(p)
            ids.append(idx)
    return ids, result


def x2_to_x(points, bondlength = 1.5):
    dmat = cdist(points, points)
    ids = cluskit.cluster._rank_fps(points, K = None, greedy =False)
    dmat = np.triu(dmat[ids, :][:, ids])
    remaining_ids = np.all((dmat > bondlength) | (dmat == 0), axis =0)
    return ids[remaining_ids == 1]

@explicit_serialize
class AdsorbateEliminationTask(FiretaskBase):
    """ 
    Task to eliminate too close adsorbate atoms.

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
                site_ids_list.append(adsorbate.get("site_ids_list", []))
            adsorbate_atoms = atoms[np.array(adsorbate_ids, dtype = int)]

            cluster_ids = []
            for nanocluster in source_simulation["nanoclusters"]:
                cluster_ids.extend(nanocluster["atom_ids"])
            cluster_atoms = atoms[np.array(cluster_ids, dtype = int)]

            adsorbate_positions = adsorbate_atoms.get_positions()
            
            #kept_positions = cluskit.utils.x2_to_x(adsorbate_positions.reshape((-1,3)), bondlength = bond_length)
            if n_remaining:
                remaining_ids = cluskit.cluster._rank_fps(adsorbate_positions, K = n_remaining, greedy =False)
            elif bond_length:
                remaining_ids = x2_to_x(adsorbate_positions, bondlength = bond_length)

            else:
                logging.warning("give either argument bond_length or n_remaining")

            #remaining_ids = x2_to_x(adsorbate_positions.reshape((-1,3)), bondlength = bond_length)
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
            #update_spec["simulations"][str(simulation_id)] = dct
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
    """ 
    Task to cover cluster fully with adsorbates based on a type of site.

    Args:
        adsorbate_name  (str) : Adsorbate atom name to be placed
                                on all sites found.
        adsite_types  (list of str) : Can be "top", "bridge" or "hollow".
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
            descriptor_setup = dscribe.descriptors.SOAP(atomic_numbers = all_atomtypes, 
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
                    dct["adsorbates"].append(dict({"atom_ids" : adsorbate_ids, "reference_id" : reference_id}))
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
            #update_spec["simulations"][str(simulation_id)] = dct
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
    """
    # scheme for coverage ladder workflow
    # pure layout of logic

        # ladder variables: l = number of minimum energies selected
        #                   d = maximum depth of branches (convergence criterion)
        #                   n0 = starting number of adsorbates
        #                   structure0 = starting structure with adsorbates

        # DAG dictionary of jobs, simulation_ids

        # dictionary of jobs per number of adsorbates (key = n, value = dict of (simulation_id : energy))

        # tracing variable how many adsorbates on the surface: n_adsorbates

        # if root calculation
        # decision to go up or down the ladder
        # based on adsorption energy

        # similarity of occupied / unoccupied sites
        # 
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
        n_adsorbates = fw_spec["temp"]["n_adsorbates"]
        l = self["l"]
        d = self["d"] # TODO implement stop criterion d
        is_new_root = fw_spec["temp"]["is_new_root"]
        n_adsorbates_root = fw_spec["temp"]["n_adsorbates_root"]
        n_adsorbates = fw_spec["temp"]["n_adsorbates"]
        open_branches = fw_spec["temp"]["open_branches"]
        root_history = fw_spec["temp"]["root_history"]
        step_history = fw_spec["temp"]["step_history"]


        lowenergy_calc_ids, energies = self.find_lowest_energy_structures(l, calc_ids, energies)

        lowest_energy, idx = np.array(energies).min(), np.argmin(np.array(energies))
        lowest_idx = lowenergy_calc_ids[idx]

        is_new_root = self.check_new_root(is_new_root, branch_dct,
                                       lowest_idx, ne_dct, n_adsorbates_root, n_adsorbates)
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
            dg_diff, assignment = self.compute_dg_diff(lowest_idx, branch_dct, ne_dct)
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
        return FWAction(update_spec=fw_spec)

    def find_lowest_energy_structures(self, l, calc_ids, energies):
        ids = np.argsort(energies)[:l] 
        calc_ids = np.array(calc_ids)[ids].tolist()
        energies = np.array(energies)[ids].tolist()
        return calc_ids, energies

    def level_check_new_root(self, idx, ne_dct, n_adsorbates_root, n_adsorbates):
        """
        Has to compare with root, if on the same coverage level.
        """

        if n_adsorbates == n_adsorbates_root:
            pass
        else:
            return False
        # requires dictionary of energies
        energies = ne_dct[n_adsorbates]
        is_new_root = energies[idx] <= min(energies.values())
        return is_new_root

    def check_new_root(self, is_new_root, branch_dct, lowest_idx, ne_dct, n_adsorbates_root, n_adsorbates):
        """
        checks if current lowest-energy configuration is eligible for new root
        """
        # first iteration after new root has been set
        if is_new_root:
            dg_diff, assignment = self.compute_dg_diff(lowest_idx, branch_dct, ne_dct)
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

    def compute_dg_diff(self, idx, branch_dct, ne_dct):
        # search where idx appears in branch_dct
        #print("IDX", idx, type(idx))
        for key, ids in zip(branch_dct.keys(), branch_dct.values()):
            if int(idx) in list(ids):
                parent_idx = key
                break
        # get energies
        # energy = find(idx, ne_dct)
        # parent_energy = find(parent_idx, ne_dct)
        
        #print("parent IDX", parent_idx, type(idx))

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
            dg_diff = parent_energy - energy
            assignment = "parent"
        else:
            dg_diff = energy - parent_energy
            assignment = "child"
        return dg_diff, assignment

    def decide_next_branch(self, open_branches, ne_dct):
        next_branch = open_branches.pop(0)
        direction = next_branch[1]
        calc_ids = next_branch[0]
        for n_adsorbates, v in ne_dct.items():
            if calc_ids[0] in v.keys():
                new_n_adsorbates = n_adsorbates
        return direction, calc_ids, new_n_adsorbates


@explicit_serialize
class GatherLadderTask(FiretaskBase):
    """
    Task for Gathering properties for coverage ladder
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


        print("COMPARISON ANALYSIS IDS, AND ORDERED")
        print(calc_ids)
        print(analysis_ids)
        print(reordered_analysis_ids)
        analysis_ids =  reordered_analysis_ids

        simulations = fetch_simulations(fw_spec["extdb_connect"], analysis_ids)
        energies = []
        for idx in analysis_ids:
            simulation = simulations[str(idx)]
            energies.append(simulation["output"]["total_energy"]) 

        # find calc_id in calc_parents and replace with analysis_id
        #print(calc_ids, len(calc_ids))
        #print(analysis_ids, len(analysis_ids))
        #print("calc_parents")
        #print(calc_parents)
        for _, children in calc_parents.items():
            for i, child in enumerate(children):
                if child in calc_ids:
                    idx = calc_ids.index(child)
                    children[i] = analysis_ids[idx]
                    #print("replacing ", calc_ids[idx], "with ", analysis_ids[idx])
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
        fw_spec["temp"]["analysis_ids"] = [] # TODO check if push from Analysis requires array
        fw_spec["temp"]["calc_parents"] = calc_parents
        fw_spec.pop("_category")
        fw_spec.pop("name")
        return FWAction(update_spec=fw_spec)


@explicit_serialize
class AddRemoveLadderTask(FiretaskBase):

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

    def add_adsorbate(self, simulation, bond_length, db, ranking_metric = "proximity", k = 7, 
        adsorbate_name = "H", reference_id = -1, workflow_id = -1):
        atoms_dict = simulation["atoms"]
        atoms = atoms_dict_to_ase(atoms_dict)
        cluster_atoms, adsorbate_atoms, site_ids_list, site_class_list, reference_ids, adsorbate_ids = self.split_nanocluster_and_adsorbates(simulation)

        # get sites of nanocluster
        cluster = cluskit.Cluster(cluster_atoms)
        adsorption_sites = cluster.get_sites(-1)

        # find "empty" sites
        ids = self.get_empty_sites(adsorption_sites, adsorbate_atoms.get_positions(), bond_length= bond_length)
        #print("get_empty_sites", ids, len(ids))

        # based on bond distance (given)
        # TODO similarity of environment
        # optional: exclude similar sites based on similarity threshold

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
        #pp(remaining_ids)

        # construct structures
        adsorbate_lst = adsorbate_pos_to_atoms_lst(adsorption_sites[remaining_ids], adsorbate_name)

        # TODO delete visualizer
        #added_adsorbates = ase.Atoms("X" * len(remaining_ids), adsorption_sites[remaining_ids])
        #view(atoms + added_adsorbates)
        ####

        # create new simulations
        # loop over each adsorbate
        simulations_chunk_list = []
        new_calc_ids = []
        for adsorbate in adsorbate_lst:
            # adsites_dict
            joint_atoms, cluster_ids, adsorbate_ids = self.join_cluster_adsorbate(atoms, adsorbate)
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
            # update_spec["simulations"][str(simulation_id)] = dct
            new_calc_ids.append(simulation_id)

        db["simulations"].insert_many(simulations_chunk_list)
        return new_calc_ids

    def remove_adsorbate(self, simulation, bond_length, db, ranking_metric = "similarity", k = 7, adsorbate_name = "H", workflow_id = -1):
        atoms_dct = simulation["atoms"]
        atoms = atoms_dict_to_ase(atoms_dict)
        cluster_atoms, adsorbate_atoms, site_ids_list, site_class_list, reference_ids, adsorbate_ids = self.split_nanocluster_and_adsorbates(
            simulation)


        # get descriptor of adsorbate atoms
        cluster = cluskit.Cluster(cluster_atoms + adsorbate_atoms)
        #cluster.get_sites(-1)
        desc = cluster.get_cluster_descriptor()[len(cluster_atoms) :]
        desc.shape


        # optional: exclude similar adsorbates based on similarity threshold

        # rank adsorbates based on
        # a) soap similarity
        # b) proximity
        if ranking_metric == "similarity":
            ids = cluskit.cluster._rank_fps(desc, K=None, greedy=False)
        else:
            ids = self.x2_to_x(adsorbate_atoms.get_positions(), bond_length = bond_length)

        # keep first k adsorbates
        if len(ids) < k:
            remaining_ids = ids
        else:
            remaining_ids = ids[:k]
        #pp(remaining_ids)

        # construct structures
        adsorbate_lst = adsorbate_pos_to_atoms_lst(adsorbate_atoms.get_positions()[remaining_ids], adsorbate_name)

        # TODO delete visualizer
        removed_adsorbates = ase.Atoms("Xe" * len(remaining_ids), adsorbate_atoms.get_positions()[remaining_ids])
        view(atoms + removed_adsorbates)
        ####

        # create new simulations
        # loop over each adsorbate
        simulations_chunk_list = []
        new_calc_ids = []


        for idx in remaining_ids:
            joint_atoms = cluster_atoms
            dct = deepcopy(simulation)
            dct["adsorbates"] = []
            for i, adsorbate_id in zip(np.arange(len(adsorbate_ids)), adsorbate_ids):
                if idx == i:
                    # removed atom is excluded
                    pass
                else:
                    adsorbate = adsorbate_atoms[idx]
                    joint_atoms, cluster_ids, adsorbate_ids = self.join_cluster_adsorbate(joint_atoms, adsorbate)
                    dct["adsorbates"].append(dict({"atom_ids": adsorbate_ids, "reference_id": reference_ids[idx],
                                                   "site_class": site_class_list[idx], "site_ids": site_ids_list[idx]}))
            # info about surface atoms not there

            # adsites_dict
            joint_atoms_dict = ase_to_atoms_dict(joint_atoms)

            # calculation originated from this:
            dct["source_id"] = simulation["_id"]
            dct["workflow_id"] = workflow_id
            dct["atoms"] = joint_atoms_dict
            dct["operations"] = [dict({"add_adsorbate": -1})]
            # empty previous input
            dct["inp"] = {}
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
            # update_spec["simulations"][str(simulation_id)] = dct
            new_calc_ids.append(simulation_id)

        db["simulations"].insert_many(simulations_chunk_list)
        return new_calc_ids

    def split_nanocluster_and_adsorbates(self, simulation):
        adsorbate_ids = []
        reference_ids = []
        site_class_list = []
        site_ids_list = []
        atoms_dict = simulation["atoms"]
        atoms = atoms_dict_to_ase(atoms_dict)

        for adsorbate in simulation["adsorbates"]:
            adsorbate_ids.extend(adsorbate["atom_ids"])
            reference_ids.append(adsorbate["reference_id"])
            site_class_list.append(adsorbate.get("site_class", ""))
            site_ids_list.append(adsorbate.get("site_ids_list", []))
        adsorbate_atoms = atoms[np.array(adsorbate_ids, dtype=int)]

        cluster_ids = []
        for nanocluster in simulation["nanoclusters"]:
            cluster_ids.extend(nanocluster["atom_ids"])
        cluster_atoms = atoms[np.array(cluster_ids, dtype=int)]

        return cluster_atoms, adsorbate_atoms, site_ids_list, site_class_list, reference_ids, adsorbate_ids

    def x2_to_x(self, points, bond_length = 1.5):
        dmat = cdist(points, points)
        ids = cluskit.cluster._rank_fps(points, K = None, greedy =False)
        dmat = np.triu(dmat[ids, :][:, ids])
        remaining_ids = np.all((dmat > bondlength) | (dmat == 0), axis =0)
        return ids[remaining_ids == 1]

    def get_empty_sites(self, site_positions, adsorbate_positions, bond_length = 1.5):
        """
        ordered by distance to next adsorbate

        :param site_positions:
        :param adsorbate_positions:
        :param bondlength:
        :return:
        """
        dmat = cdist(site_positions, adsorbate_positions)
        mindist = dmat.min(axis = 1)
        ids = np.argsort(-mindist)
        mindist = mindist[ids]
        remaining_ids = mindist > bond_length
        return ids[remaining_ids == 1]

    def join_cluster_adsorbate(self, cluster, adsorbate):
        joint_atoms = cluster + adsorbate
        cluster_ids = list(range(len(cluster)))
        adsorbate_ids = list(range(len(cluster_ids), len(joint_atoms)))

        return joint_atoms, cluster_ids, adsorbate_ids

@explicit_serialize
class NewLadderRootTask(FiretaskBase):
    """
    Task to start CoverageLadder Workflow
    """

    _fw_name = 'NewLadderRootTask'
    required_params = ['start_ids',]
    optional_params = ['initial_direction']

    def run_task(self, fw_spec):
        start_ids = self["start_ids"]
        initial_direction = self["initial_direction"]
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
        fw_spec.pop("_category")
        fw_spec.pop("name")
        return FWAction(update_spec=fw_spec)


def start_coverage_ladder(start_ids, initial_direction = 1):
    firetask1  = NewLadderRootTask(start_ids = start_ids, initial_direction = initial_direction)
    fw = Firework([firetask1], spec = {'_category' : "small", 'name' : 'NewLadderRootTask'},
             name = 'NewLadderRootWork')
    return fw

def add_remove_adsorbate(bond_length, k = 7, ranking_metric = "similarity"):
    firetask1  = AddRemoveLadderTask(bond_length = bond_length, k = k, ranking_metric = ranking_metric)
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'AddRemoveLadderTask'},
             name = 'AddRemoveLadderWork')    
    return fw

def gather_ladder():
    firetask1  = GatherLadderTask()
    fw = Firework([firetask1], spec = {'_category' : "small", 'name' : 'GatherLadderTask'},
             name = 'GatherLadderWork')    
    return fw

def step_coverage_ladder(l = 2, d = 4):
    firetask1  = CoverageLadderTask(l = l, d = d)
    fw = Firework([firetask1], spec = {'_category' : "small", 'name' : 'CoverageLadderTask'},
             name = 'CoverageLadderWork')    
    return fw


def get_per_type_coverage(reference_energy = 0.0, adsorbate_name='H', adsite_types = ["top", "bridge", "hollow"], 
        descriptor = "soap", descriptor_params = {"nmax" : 9, "lmax" :6, "rcut" : 5.0, 
            "crossover" : True, "sparse" : False}):
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
    firetask1  = AdsorbateEliminationTask(adsorbate_name = adsorbate_name, bond_length = bond_length)
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'AdsorbateEliminationTask'},
             name = 'AdsorbateEliminationWork')
    return fw

def eliminate_closest(adsorbate_name, n_remaining):
    firetask1  = AdsorbateEliminationTask(adsorbate_name = adsorbate_name, n_remaining = n_remaining)
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'AdsorbateEliminationTask'},
             name = 'AdsorbateEliminationWork')
    return fw
