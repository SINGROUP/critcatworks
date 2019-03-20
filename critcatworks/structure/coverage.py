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

from critcatworks.database import atoms_dict_to_ase, ase_to_atoms_dict
from critcatworks.database import read_descmatrix, write_descmatrix
from critcatworks.database.extdb import update_simulations_collection
from critcatworks.database.extdb import fetch_simulations
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
            if bond_length:
                remaining_ids = x2_to_x(adsorbate_positions, bondlength = bond_length)
            elif n_remaining:
                remaining_ids = cluskit.cluster._rank_fps(points, K = n_remaining, greedy =False)

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
            
        #update_spec["temp"]["descmatrix"] = descmatrix
        # saves descmatrix as a path to a numpy array
        update_spec["temp"]["descmatrix"] = write_descmatrix(descmatrix)
        update_spec["temp"]["calc_ids"] = new_calc_ids
        print("spec bytes")
        print(sys.getsizeof(update_spec))
        #print("descmatrix bytes")
        #print(descmatrix.nbytes)
        #print(descmatrix.shape)
        import json
        #print(len(json.dumps(update_spec["temp"]["descmatrix"])))
        #print(len(json.dumps(update_spec["simulations"])))
        print("lenght of spec in json", len(json.dumps(update_spec)))
        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)


@explicit_serialize
class CoverageLadderTask(FiretaskBase):
    """ 
    Task to go up or down the coverage ladder.

    Args:

    """

    _fw_name = 'CoverageLadderTask'
    required_params = []
    optional_params = []

    def run_task(self, fw_spec):

        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec) 


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
