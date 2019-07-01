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
from critcatworks.database import adsorbate_pos_to_atoms_lst
from critcatworks.database import join_cluster_adsorbate
from critcatworks.database.extdb import gather_all_atom_types
from critcatworks.database.extdb import update_simulations_collection
from critcatworks.database.extdb import get_external_database, _query_id_counter_and_increment
from critcatworks.database.extdb import fetch_simulations


@explicit_serialize
class AdsiteCreationTask(FiretaskBase):
    """ 
    Firetask to determine adsorption site structures. The adsorbate
    can only be a single atom. For a molecular adsorbate (binding
    to only one site), please use the Firetask 
    MonodentateAdsiteCreationTask.
    On a set of nanoclusters, a combination of top/bridge/hollow sites
    is classified and populated. Hence, the simulation collection
    of the mongodb database is updated with structures which each
    contain a single adsorbate on a nanocluster.

    Besides, a descriptor matrix is written to a file for later use.

    Args:
        reference_energy (float) :  reference energy for the adsorbate. Can be the
                                    total energy of the isolated adsorbate molecule
                                    or a different reference point
        adsorbate_name (str) :  element symbold of the adsorbed atom
        adsite_types (list) :   adsorption site types, can contain any combination of
                                "top", "bridge", "hollow"
                
    Returns:
        FWAction : Firework action, updates fw_spec
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
            # TODO allow for other descriptors than default
            # TODO read descriptor_params
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
                    simulations_chunk_list.append(dct)

                    # update internal workflow data
                    simulation_id = dct["_id"]
                    new_calc_ids.append(simulation_id)

            db["simulations"].insert_many(simulations_chunk_list)
        
        descmatrix = np.array(desc_lst)

        # saves descmatrix as a path to a numpy array
        update_spec["temp"]["descmatrix"] = write_descmatrix(descmatrix)
        update_spec["temp"]["calc_ids"] = new_calc_ids

        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)

@explicit_serialize
class MonodentateAdsiteCreationTask(FiretaskBase):
    """ 
    Firetask to determine adsorption site structures. The adsorbate
    is specified as a molecule with an anchor X which binds to
    only one site. The position of the binding atom is exactly
    determined by the adsorption vector of the adsorption site
    and the bond length between X and the closest, later binding,
    atom. On a set of nanoclusters, a combination of 
    top/bridge/hollow sites is classified and populated. 
    Hence, the simulation collection of the mongodb database 
    is updated with structures which each
    contain a single adsorbate on a nanocluster.

    Besides, a descriptor matrix is written to a file for later use.

    Args:
        reference_energy (float) :  reference energy for the adsorbate. Can be the
                                    total energy of the isolated adsorbate molecule
                                    or a different reference point
        adsorbate (dict) :  adsorbed molecule as atoms dict. Contains an "X" dummy atom
                            which indicates the anchor point to the nanocluster
        adsite_types (list) :   adsorption site types, can contain any combination of
                                "top", "bridge", "hollow"
                
    Returns:
        FWAction : Firework action, updates fw_spec
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
            # TODO allow for other descriptors than default
            # TODO read descriptor_params
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
                    simulations_chunk_list.append(dct)

                    # update internal workflow data
                    simulation_id = dct["_id"]
                    new_calc_ids.append(simulation_id)
        
            db["simulations"].insert_many(simulations_chunk_list)
        descmatrix = np.array(desc_lst)

        # saves descmatrix as a path to a numpy array
        update_spec["temp"]["descmatrix"] = write_descmatrix(descmatrix)
        update_spec["temp"]["calc_ids"] = new_calc_ids

        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)


@explicit_serialize
class MonodentateUniqueAdsiteCreationTask(FiretaskBase):
    """ 
    Firetask to determine unique adsorption site structures. The 
    adsorbate is specified as a molecule with an anchor X which 
    binds to only one site. The position of the binding atom is 
    exactly determined by the adsorption vector of the adsorption 
    site and the bond length between X and the closest, later 
    binding, atom. On a set of nanoclusters, a combination of 
    top/bridge/hollow sites is classified.
    The sites are ranked for their dissimilarity based on a 
    local descriptor (default soap), then uniqueness is defined
    by a threshold parameter. Only the unique sites are kept and
    populated with the adsorbate. 

    Hence, the simulation collection of the mongodb database 
    is updated with structures which each contain a single adsorbate 
    on a nanocluster.

    Besides, a descriptor matrix is written to a file for later use.

    Args:
        reference_energy (float) :  reference energy for the adsorbate. Can be the
                                    total energy of the isolated adsorbate molecule
                                    or a different reference point
        adsorbate (dict) :  adsorbed molecule as atoms dict. Contains an "X" dummy atom
                            which indicates the anchor point to the nanocluster
        adsite_types (list) :   adsorption site types, can contain any combination of
                                "top", "bridge", "hollow"
        threshold (float) :     threshold of similarity metric between the local structures
                                of the adsorption sites. Only sites which are more dissimilar 
                                than the given threshold are computed

    Returns:
        FWAction : Firework action, updates fw_spec
    """
    _fw_name = 'MonodentateUniqueAdsiteCreationTask'
    required_params = ['reference_energy', 'adsorbate', 'adsite_types', 'threshold']
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
            # TODO allow for other descriptors than default
            # TODO read descriptor_params
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

                # reduce to unique sites per site type
                unique_ids = cluster.get_unique_sites(sitetype = adsite_type_int)
                adsorbate_lst = [adsorbate_lst[idx] for idx in unique_ids]

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
                    simulations_chunk_list.append(dct)

                    # update internal workflow data
                    simulation_id = dct["_id"]
                    new_calc_ids.append(simulation_id)
        
            db["simulations"].insert_many(simulations_chunk_list)
        descmatrix = np.array(desc_lst)

        # saves descmatrix as a path to a numpy array
        update_spec["temp"]["descmatrix"] = write_descmatrix(descmatrix)
        update_spec["temp"]["calc_ids"] = new_calc_ids

        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)



@explicit_serialize
class AdsiteRankTask(FiretaskBase):
    """ 
    Firetask to determine ranking of adsorption site structures.
    Based on a previously calculated descriptor matrix (in
    the fw_spec as a path to a file), the adsorbate structure 
    sites are ranked are ranked for their dissimilarity.

    Hence, the calc_ids, the structures to be simulated are
    reordered. 

    Besides, the descriptor matrix is reordered and written 
    to a file for later use.

    Args:

    Returns:
        FWAction : Firework action, updates fw_spec
    """

    _fw_name = 'AdsiteRankTask'
    required_params = []
    optional_params = []

    def run_task(self, fw_spec):
        logging.debug(fw_spec)
        calc_ids = fw_spec["temp"]["calc_ids"]
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
    """ 
    Firetask to determine adsorption site structures. The adsorbate
    can only be a single atom. For a molecular adsorbate (binding
    to only one site), please use the Firetask 
    MonodentateAdsiteCreationTask.
    On a set of nanoclusters, a combination of top/bridge/hollow sites
    is classified and populated. Hence, the simulation collection
    of the mongodb database is updated with structures which each
    contain a single adsorbate on a nanocluster.

    Besides, a descriptor matrix is written to a file for later use.

    Args:
        reference_energy (float) :  reference energy for the adsorbate. Can be the
                                    total energy of the isolated adsorbate molecule
                                    or a different reference point
        adsorbate_name (str) :  element symbold of the adsorbed atom
        adsite_types (list) :   adsorption site types, can contain any combination of
                                "top", "bridge", "hollow"
        descriptor (str) :  type of descriptor to be used. For a list of
                            descriptors, see the documentation of dscribe
                            Defaults to 'soap'
        descriptor_params (dict) : descriptor parameters
                    
    Returns:
        FWAction : Firework action, updates fw_spec
    """
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
    """ 
    Firetask to determine adsorption site structures. The adsorbate
    is specified as a molecule with an anchor X which binds to
    only one site. The position of the binding atom is exactly
    determined by the adsorption vector of the adsorption site
    and the bond length between X and the closest, later binding,
    atom. On a set of nanoclusters, a combination of 
    top/bridge/hollow sites is classified and populated. 
    Hence, the simulation collection of the mongodb database 
    is updated with structures which each
    contain a single adsorbate on a nanocluster.

    Besides, a descriptor matrix is written to a file for later use.

    Args:
        reference_energy (float) :  reference energy for the adsorbate. Can be the
                                    total energy of the isolated adsorbate molecule
                                    or a different reference point
        adsorbate (dict) :  adsorbed molecule as atoms dict. Contains an "X" dummy atom
                            which indicates the anchor point to the nanocluster
        adsite_types (list) :   adsorption site types, can contain any combination of
                                "top", "bridge", "hollow"
        descriptor (str) :  type of descriptor to be used. For a list of
                            descriptors, see the documentation of dscribe
                            Defaults to 'soap'
        descriptor_params (dict) : descriptor parameters

    Returns:
        Firework : Firework MonodentateAdsiteCreationWork
    """
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

def get_monodentate_unique_adsites(reference_energy = 0.0, adsorbate = {}, adsite_types = ["top", "bridge", "hollow"],
        threshold = 0.001, descriptor = "soap", 
        descriptor_params = {"nmax" : 9, "lmax" :6, "rcut" : 5.0, 
            "crossover" : True, "sparse" : False}):
    """ 
    Firework to determine unique adsorption site structures. The 
    adsorbate is specified as a molecule with an anchor X which 
    binds to only one site. The position of the binding atom is 
    exactly determined by the adsorption vector of the adsorption 
    site and the bond length between X and the closest, later 
    binding, atom. On a set of nanoclusters, a combination of 
    top/bridge/hollow sites is classified.
    The sites are ranked for their dissimilarity based on a 
    local descriptor (default soap), then uniqueness is defined
    by a threshold parameter. Only the unique sites are kept and
    populated with the adsorbate. 

    Hence, the simulation collection of the mongodb database 
    is updated with structures which each contain a single adsorbate 
    on a nanocluster.

    Besides, a descriptor matrix is written to a file for later use.

    Args:
        reference_energy (float) :  reference energy for the adsorbate. Can be the
                                    total energy of the isolated adsorbate molecule
                                    or a different reference point
        adsorbate (dict) :  adsorbed molecule as atoms dict. Contains an "X" dummy atom
                            which indicates the anchor point to the nanocluster
        adsite_types (list) :   adsorption site types, can contain any combination of
                                "top", "bridge", "hollow"
        threshold (float) :     threshold of similarity metric between the local structures
                                of the adsorption sites. Only sites which are more dissimilar 
                                than the given threshold are computed
        descriptor (str) :  type of descriptor to be used. For a list of
                            descriptors, see the documentation of dscribe
                            Defaults to 'soap'
        descriptor_params (dict) : descriptor parameters

    Returns:
        Firework : Firework MonodentateUniqueAdsiteCreationWork
    """
    firetask1  = MonodentateUniqueAdsiteCreationTask(
        adsorbate = adsorbate, 
        adsite_types = adsite_types,
        reference_energy = reference_energy,
        descriptor = descriptor,
        descriptor_params = descriptor_params,
        threshold = threshold,
        )
    fw = Firework([firetask1], 
        spec={'_category' : "medium", 'name' : 'MonodentateUniqueAdsiteCreationTask'},
        name = 'MonodentateUniqueAdsiteCreationWork'
        )
    return fw

def rank_adsites():
    """    
    Firework to determine ranking of adsorption site structures.
    Based on a previously calculated descriptor matrix (in
    the fw_spec as a path to a file), the adsorbate structure 
    sites are ranked are ranked for their dissimilarity.

    Hence, the calc_ids, the structures to be simulated are
    reordered. 

    Besides, the descriptor matrix is reordered and written 
    to a file for later use.

    Args:

    Returns:
        Firework : Firework AdsiteRankWork
    """
    firetask1  = AdsiteRankTask()
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'AdsiteRankTask'},
             name = 'AdsiteRankWork')
    return fw
