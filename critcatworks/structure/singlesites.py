from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time, pathlib, sys
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import ase, ase.io
import clusgeo
import dscribe
import numpy as np
import logging

from critcatworks.database import atoms_dict_to_ase, ase_to_atoms_dict

def adsorbate_pos_to_atoms_dict(structure, adspos, adsite_type):
    """
    out of use. TODO delete.
    """
    atoms_lst = []
    ads_structures_dict = []
    for adatom, idx in zip(adspos, range(len(adspos))):
        logging.debug(adsite_type)
        logging.debug(adatom)
        logging.debug(adatom.shape)
        adatom = ase.Atoms(symbols=adsite_type, positions=adatom.reshape((1,3)))
        clus_ads = structure + adatom
        atoms_lst.append(clus_ads)
        clus_ads_dict = ase_to_atoms_dict(clus_ads)
        ads_structures_dict.append(clus_ads_dict)

    return ads_structures_dict

def adsorbate_pos_to_atoms_lst(adspos, adsite_type):
    """
    works with only one adsorbate atom.
    TODO generalize cluskit to return list of adsorbates
    already in ase format.
    Makes this function superfluous.
    """
    atoms_lst = []
    ads_structures_dict = []
    for adsorbate, idx in zip(adspos, range(len(adspos))):
        logging.debug(adsite_type)
        logging.debug(adatom)
        logging.debug(adatom.shape)
        adatom = ase.Atoms(symbols=adsite_type, positions=adatom.reshape((1,3)))
        atoms_lst.append(clus_ads)
    return atoms_lst

def join_cluster_adsorbate(cluster, adsorbate):
    joint_atoms = cluster + adsorbate
    cluster_ids = list(range(len(cluster)))
    adsorbate_ids = list(range(len(cluster_ids), len(joint_atoms)))

    return joint_atoms, cluster_ids, adsorbate_ids

def gather_all_atom_types(calc_ids, simulations_collection):
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
        simulations = fw_spec["simulations"]
        workflow_id = fw_spec.get("workflow", {"_id" : -1 }).get("_id", -1)
        update_spec = fw_spec

        logging.debug(fw_spec)
        desc_lst = []
        new_calc_ids = []

        # create reference of adsorbate in order to store its total energy
        # for later constructing adsorption energies
        reference_simulation = update_simulations_collection(atoms = {}, 
            source_id = -1, workflow_id = workflow_id, 
            nanoclusters = [], adsorbates = [], substrates = [], 
            operations = [""], inp = {"adsorbate_name" : adsorbate_name}, 
            output = {"total_energy" : reference_energy},)
        reference_id = reference_simulation["_id"]


        all_atomtypes = gather_all_atom_types(calc_ids, simulations)

        # looping over nc atoms 
        for idx, atoms in enumerate(calc_ids):
            logging.debug(atoms)

            ##
            # get source simulation
            source_simulation = simulations[str(calc_id)]
            # entries for adsorbate dictionary
            ads_pos_lst = []

            # running clusgeo on cluster
            cluster = clusgeo.ClusGeo(atoms)
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


                adsorbate_lst = adsorbate_pos_to_atoms_lst(adspos, adsite_type)

                #loop over each adsorbate
                for adsorbate, surface_atoms in zip(adsorbate_lst, sites_surface_atoms):

                    #adsites_dict 
                    joint_atoms, cluster_ids, adsorbate_ids = join_cluster_adsorbate(cluster, adsorbate)
                    joint_atoms_dict = ase_to_atoms_dict(joint_atoms)

                    # update external database
                    dct = source_simulation.copy()
                    # calculation originated from this:
                    dct["source_id"] = calc_id
                    dct["workflow_id"] = workflow_id
                    dct["atoms"] = joint_atoms_dict
                    dct["operations"] = [{"add_adsorbate" : 1}]
                    dct["adsorbates"].append({"atom_ids" : adsorbate_ids, "reference_id" : reference_id})
                    dct["inp"]["adsite_type"] = adsite_type
                    dct["inp"]["adsorbate_name"] = adsorbate_name
                    dct["output"]["surface_atoms"] = surface_atoms

                    logging.info(dct)
                    simulation = update_simulations_collection(dct)

                    # update internal workflow data
                    simulation_id = simulation["_id"]
                    update_spec["simulations"][str(simulation_id)] = dct
                    new_calc_ids.append(simulation_id)



            # move coverage to other task


                ### full coverage addition : make coverage structure anyway and add to simulation
                # get structure with all adsorbates
                #ads_pos_lst.append(adsites)

            # get structure with all adsorbates, different sites combined
            #clus_cov = atoms.copy()
            #for pos in ads_pos_lst:
            #    adatoms = ase.Atoms(symbols=[adsorbate_name] * pos.shape[0], positions=pos)
            #    clus_cov += adatoms
            #clus_cov_dict = ase_to_atoms_dict(clus_cov)
            #coverage_structures.append(clus_cov_dict)
            #coverage_id_dict[str(idx)] = []
        
        descmatrix = np.array(desc_lst)

            
        update_spec["temp"]["descmatrix"] = descmatrix
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
        descmatrix = fw_spec["temp"]["descmatrix"]
        descmatrix = np.array(descmatrix)
        logging.info("DESCRIPTOR matrix attributes")
        logging.info(descmatrix.shape)
        logging.info(np.sum(descmatrix))
        fps_ranking = clusgeo.cluster._rank_fps(descmatrix, K = None, greedy =False, is_safe = True)
        update_spec = fw_spec
        update_spec["temp"]["fps_ranking"] = fps_ranking
        update_spec.pop("_category")
        update_spec.pop("name")

        return FWAction(update_spec=update_spec)






def get_adsites(adsorbate_name='H', adsite_types = ["top", "bridge", "hollow"]):
    firetask1  = AdsiteCreationTask(
        adsorbate_name=adsorbate_name, 
        adsite_types = adsite_types,
        )
    fw = Firework([firetask1], 
        spec={'adsorbate_name' : adsorbate_name, 
        'adsite_types' : adsite_types,
        '_category' : "lightweight",
        'name' : 'AdsiteCreationTask'},
        name = 'AdsiteCreationWork'
        )
    return fw


def rank_adsites():
    firetask1  = AdsiteRankTask()
    fw = Firework([firetask1], spec = {'_category' : "medium", 'name' : 'AdsiteRankTask'},
             name = 'AdsiteRankWork')
    return fw
