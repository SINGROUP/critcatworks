from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time, pathlib, sys, copy
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import ase, ase.io
import cluskit
import dscribe
import numpy as np
import logging

from critcatworks.database import atoms_dict_to_ase, ase_to_atoms_dict
from critcatworks.database import read_descmatrix, write_descmatrix
from critcatworks.database.extdb import update_simulations_collection
from critcatworks.database.extdb import get_external_database, _query_id_counter_and_increment
from critcatworks.database.extdb import fetch_simulations

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

def join_cluster_adsorbate(cluster, adsorbate):
    joint_atoms = cluster + adsorbate
    cluster_ids = list(range(len(cluster)))
    adsorbate_ids = list(range(len(cluster_ids), len(joint_atoms)))

    return joint_atoms, cluster_ids, adsorbate_ids

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
class AdsiteCreationTask(FiretaskBase):
    """ 
    Task to determine adsorption site structures.

    Args:
        adsorbate_name  (str) : Adsorbate atom name to be placed
                                on all sites found.
        adsite_types  (list of str) : Can be "top", "bridge" or "hollow".
    """

    _fw_name = 'AdsiteCreationTask'
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
        db = get_external_database(fw_spec["extdb_connect"])

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
            simulations_chunk_list = []
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

                    #adsites_dict 
                    joint_atoms, cluster_ids, adsorbate_ids = join_cluster_adsorbate(atoms, adsorbate)
                    joint_atoms_dict = ase_to_atoms_dict(joint_atoms)

                    # update external database
                    dct = copy.deepcopy(source_simulation)
                    # calculation originated from this:
                    dct["source_id"] = calc_id
                    dct["workflow_id"] = workflow_id
                    dct["atoms"] = joint_atoms_dict
                    dct["operations"] = [dict({"add_adsorbate" : 1})]
                    dct["adsorbates"].append(dict({"atom_ids" : adsorbate_ids, "reference_id" : reference_id}))
                    # empty previous input
                    dct["inp"] = {}
                    dct["inp"]["adsite_type"] = adsite_type
                    dct["inp"]["adsorbate_name"] = adsorbate_name
                    # empty previous output
                    dct["output"] = {}
                    dct["output"]["surface_atoms"] = surface_atoms.tolist()

                    # getting only id for uploading simulations in chunks
                    dct["_id"] = _query_id_counter_and_increment('simulations', db)
                    #simulation = update_simulations_collection(extdb_connect = fw_spec["extdb_connect"], **dct)
                    simulations_chunk_list.append(dct)

                    # update internal workflow data
                    simulation_id = dct["_id"]
                    #update_spec["simulations"][str(simulation_id)] = dct
                    new_calc_ids.append(simulation_id)

            db["simulations"].insert_many(simulations_chunk_list)
        
        descmatrix = np.array(desc_lst)

            
        #update_spec["temp"]["descmatrix"] = descmatrix
        # saves descmatrix as a path to a numpy array
        update_spec["temp"]["descmatrix"] = write_descmatrix(descmatrix)
        update_spec["temp"]["calc_ids"] = new_calc_ids

        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)

@explicit_serialize
class MonodentateAdsiteCreationTask(FiretaskBase):
    """ 
    Task to determine adsorption site structures.

    Args:
        adsorbate  (str) : Adsorbate molecule with anchor atom x in json format.

        adsite_types  (list of str) : Can be "top", "bridge" or "hollow".
    """

    _fw_name = 'MonodentateAdsiteCreationTask'
    required_params = ['reference_energy', 'adsorbate', 'adsite_types']
    optional_params = []

    def run_task(self, fw_spec):
        adsorbate_dict = self["adsorbate"]
        adsite_types = self["adsite_types"]
        reference_energy = self["reference_energy"]
        calc_ids = fw_spec["temp"]["calc_ids"]
        simulations = fetch_simulations(fw_spec["extdb_connect"], calc_ids)
        workflow_id = fw_spec.get("workflow", {"_id" : -1 }).get("_id", -1)
        update_spec = fw_spec

        logging.debug(fw_spec)
        desc_lst = []
        new_calc_ids = []
        db = get_external_database(fw_spec["extdb_connect"])

        # adsorbate atom with anchor x
        adsorbate_x = atoms_dict_to_ase(adsorbate_dict)


        # create reference of adsorbate in order to store its total energy
        # for later constructing adsorption energies
        reference_simulation = update_simulations_collection(extdb_connect = fw_spec["extdb_connect"], atoms = adsorbate_dict, 
            source_id = -1, workflow_id = workflow_id, 
            nanoclusters = [], adsorbates = [], substrates = [], 
            operations = [""], inp = {"adsorbate" : "monodentate_molecule"}, 
            output = {"total_energy" : reference_energy},)
        reference_id = reference_simulation["_id"]


        all_atomtypes = gather_all_atom_types(calc_ids, simulations)

        # looping over nc atoms 
        for idx, calc_id in enumerate(calc_ids):
            simulations_chunk_list = []
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
                _ = cluster.get_sites(adsite_type_int, distance = 0.0)
                sites_surface_atoms = cluster.site_surface_atom_ids[adsite_type_int]

                # get descriptor
                desc = cluster.get_sites_descriptor(adsite_type_int)
                for i in range(desc.shape[0]):
                    desc_lst.append(desc[i])

                adsorbate_lst = cluster.place_adsorbates(adsorbate_x, sitetype = adsite_type_int)

                #loop over each adsorbate
                for adsorbate, surface_atoms in zip(adsorbate_lst, sites_surface_atoms):

                    #adsites_dict 
                    joint_atoms, cluster_ids, adsorbate_ids = join_cluster_adsorbate(atoms, adsorbate)
                    joint_atoms_dict = ase_to_atoms_dict(joint_atoms)

                    # update external database
                    dct = copy.deepcopy(source_simulation)
                    # calculation originated from this:
                    dct["source_id"] = calc_id
                    dct["workflow_id"] = workflow_id
                    dct["atoms"] = joint_atoms_dict
                    dct["operations"] = [dict({"add_adsorbate" : 1})]
                    dct["adsorbates"].append(dict({"atom_ids" : adsorbate_ids, "reference_id" : reference_id}))
                    # empty previous input
                    dct["inp"] = {}
                    dct["inp"]["adsite_type"] = adsite_type
                    dct["inp"]["adsorbate"] = adsorbate_dict
                    # empty previous output
                    dct["output"] = {}
                    dct["output"]["surface_atoms"] = surface_atoms.tolist()

                    # getting only id for uploading simulations in chunks
                    dct["_id"] = _query_id_counter_and_increment('simulations', db)
                    #simulation = update_simulations_collection(extdb_connect = fw_spec["extdb_connect"], **dct)
                    simulations_chunk_list.append(dct)

                    # update internal workflow data
                    simulation_id = dct["_id"]
                    #update_spec["simulations"][str(simulation_id)] = dct
                    new_calc_ids.append(simulation_id)
        
            db["simulations"].insert_many(simulations_chunk_list)
        descmatrix = np.array(desc_lst)

            
        #update_spec["temp"]["descmatrix"] = descmatrix
        # saves descmatrix as a path to a numpy array
        update_spec["temp"]["descmatrix"] = write_descmatrix(descmatrix)
        update_spec["temp"]["calc_ids"] = new_calc_ids

        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)



@explicit_serialize
class AdsiteRankTask(FiretaskBase):
    """ 
    Task to determine ranking of adsorption site structures.

    """

    _fw_name = 'AdsiteRankTask'
    required_params = []
    optional_params = []

    def run_task(self, fw_spec):
        logging.debug(fw_spec)
        calc_ids = fw_spec["temp"]["calc_ids"]
        #descmatrix = fw_spec["temp"]["descmatrix"]
        #descmatrix = np.array(descmatrix)
        descmatrix = read_descmatrix(fw_spec)
        logging.info("DESCRIPTOR matrix attributes")
        logging.info(descmatrix.shape)
        logging.info(np.sum(descmatrix))
        fps_ranking = cluskit.cluster._rank_fps(descmatrix, K = None, greedy =False)

        reordered_calc_ids = np.array(calc_ids)[fps_ranking]
        reordered_descmatrix = descmatrix[fps_ranking]
        update_spec = fw_spec
        update_spec["temp"]["fps_ranking"] = fps_ranking.tolist()
        update_spec["temp"]["calc_ids"] = reordered_calc_ids.tolist()
        update_spec["temp"]["descmatrix"] = write_descmatrix(reordered_descmatrix)
        update_spec.pop("_category")
        update_spec.pop("name")

        return FWAction(update_spec=update_spec)






def get_adsites(reference_energy = 0.0, adsorbate_name='H', adsite_types = ["top", "bridge", "hollow"], 
        descriptor = "soap", descriptor_params = {"nmax" : 9, "lmax" :6, "rcut" : 5.0, 
            "crossover" : True, "sparse" : False}):
    firetask1  = AdsiteCreationTask(
        adsorbate_name=adsorbate_name, 
        adsite_types = adsite_types,
        reference_energy = reference_energy,
        descriptor = descriptor,
        descriptor_params = descriptor_params
        )
    fw = Firework([firetask1], 
        spec={'_category' : "medium", 'name' : 'AdsiteCreationTask'},
        name = 'AdsiteCreationWork'
        )
    return fw

def get_monodentate_adsites(reference_energy = 0.0, adsorbate = {}, adsite_types = ["top", "bridge", "hollow"], 
        descriptor = "soap", descriptor_params = {"nmax" : 9, "lmax" :6, "rcut" : 5.0, 
            "crossover" : True, "sparse" : False}):
    firetask1  = MonodentateAdsiteCreationTask(
        adsorbate = adsorbate, 
        adsite_types = adsite_types,
        reference_energy = reference_energy,
        descriptor = descriptor,
        descriptor_params = descriptor_params
        )
    fw = Firework([firetask1], 
        spec={'_category' : "medium", 'name' : 'MonodentateAdsiteCreationTask'},
        name = 'MonodentateAdsiteCreationWork'
        )
    return fw

def rank_adsites():
    firetask1  = AdsiteRankTask()
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'AdsiteRankTask'},
             name = 'AdsiteRankWork')
    return fw
