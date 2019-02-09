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

from critcatworks.database import atoms_dict_to_ase, ase_to_atoms_dict
from critcatworks.database.extdb import update_simulations_collection
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

@explicit_serialize
class AdsorbateEliminationTask(FiretaskBase):
    """ 
    Task to eliminate too close adsorbate atoms.

    """

    _fw_name = 'AdsorbateEliminationTask'
    required_params = ['adsorbate_name', 'bond_length']
    optional_params = []

    def run_task(self, fw_spec):
        logging.debug(fw_spec)
        adsorbate_name = self['adsorbate_name']
        bond_length = self['bond_length']

        calc_ids = fw_spec["temp"]["calc_ids"]
        new_calc_ids = []
        update_spec = fw_spec

        # divide adsorbate and nanocluster
        for idx, calc_id in enumerate(calc_ids):
            source_simulation = fw_spec["simulations"][str(calc_id)]
            atoms = atoms_dict_to_ase(source_simulation["atoms"])
            dct = copy.deepcopy(source_simulation)
            dct["source_id"] = calc_id
            dct["inp"] = {}
            dct["inp"]["bond_length"] = bond_length
            dct["inp"]["adsorbate_name"] = adsorbate_name
            dct["adsorbates"] = []
            dct["output"] = {}

            adsorbate_ids = []

            for adsorbate in source_simulation["adsorbates"]:
                adsorbate_ids.append(adsorbate["atom_ids"])
            adsorbate_atoms = atoms[np.array(adsorbate_ids, dtype = int)]

            cluster_ids = []
            for nanocluster in source_simulation["nanoclusters"]:
                cluster_ids.append(nanocluster["atom_ids"])
            cluster_atoms = atoms[np.array(cluster_ids, dtype = int)]

            adsorbate_positions = adsorbate_atoms.get_positions()
            
            kept_positions = cluskit.utils.x2_to_x(adsorbate_positions.reshape((-1,3)), bondlength = bond_length)
            n_kept = kept_positions.shape[0]
            kept_adsorbate_atoms = ase.Atoms(symbols=[adsorbate_name] * n_kept, positions=kept_positions)


            new_atoms = cluster_atoms
            for adsorbate in kept_adsorbate_atoms:
                new_atoms, cluster_ids, adsorbate_ids = join_cluster_adsorbate(atoms, adsorbate)

                #add_adsorbate += 1
                atoms, _ , adsorbate_ids = join_cluster_adsorbate(atoms, adsorbate)
                # TODO reference ID through keeping track of which atoms are removed
                dct["adsorbates"].append(dict({"atom_ids" : adsorbate_ids, "reference_id" : "NONE"}))
                # info about surface atoms not there
                # elimination function does not keep track of atom indices
            n_adsorbates = len(dct["adsorbates"]) - len(adsorbate_atoms)            
            dct["operations"] = [dict({"add_adsorbate" : n_adsorbates})]
            if n_adsorbates == 0:
                defuse_workflow = True
            else:
                defuse_workflow = False
            dct["atoms"] = ase_to_atoms_dict(new_atoms)

            # update simulation dict
            # refresh adsorbates list
            # create new simulation
            simulation = update_simulations_collection(extdb_connect = fw_spec["extdb_connect"], **dct)
            logging.info("simulation after adding single adsorbates")
            logging.info(simulation)

            # update internal workflow data
            simulation_id = simulation["_id"]
            update_spec["simulations"][str(simulation_id)] = dct
            new_calc_ids.append(simulation_id)

        update_spec["temp"]["calc_ids"] = new_calc_ids
        update_spec.pop("_category")
        update_spec.pop("name")

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
        simulations = fw_spec["simulations"]
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
            update_spec["simulations"][str(simulation_id)] = dct
            new_calc_ids.append(simulation_id)
        
        descmatrix = np.array(desc_lst)
            
        update_spec["temp"]["descmatrix"] = descmatrix
        update_spec["temp"]["calc_ids"] = new_calc_ids

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
