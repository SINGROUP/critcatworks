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
import numpy as np
import logging

from critcatworks.database import atoms_dict_to_ase

def adsorbate_pos_to_atoms_dict(structure, adspos, adsite_type):
    atoms_lst = []
    ads_structures_dict = []
    for adatom, idx in zip(adspos, range(len(adspos))):
        logging.debug(adsite_type)
        logging.debug(adatom)
        logging.debug(adatom.shape)
        adatom = ase.Atoms(symbols=adsite_type, positions=adatom.reshape((1,3)))
        clus_ads = structure + adatom
        atoms_lst.append(clus_ads)
        clus_ads_dict = clus_ads.__dict__
        ads_structures_dict.append(clus_ads_dict)

    return ads_structures_dict


def adsorbate_pos_to_descriptor(structure, adspos, adsite_type, all_atomtypes, rcut = 5.0):

    all_atomtypes = [int(i) for i in all_atomtypes]
    desc = clusgeo.environment.get_soap_sites(structure, adspos, rCut=rcut, NradBas=9, Lmax=6, 
        crossOver=True, all_atomtypes=all_atomtypes)
    return desc


@explicit_serialize
class AdsiteCreationTask(FiretaskBase):
    """ 
    Task to determine adsorption site structures.

    Args:
        reference_energy (float) :  Total energy of the free adsorbate 
                                    as a reference energy to 
                                    determine the adsorption energy 
                                    as in
                                    dE = dE(all) - dE(nc) - dE(adsorbate)
        adsorbate_name  (str) : Adsorbate atom name to be placed
                                on all sites found.
        adsite_types  (list of str) : Can be "top", "bridge" or "hollow".
    """

    _fw_name = 'AdsiteCreationTask'
    required_params = ['reference_energy', 'adsorbate_name', 'adsite_types']
    optional_params = []

    def run_task(self, fw_spec):
        reference_energy = self["reference_energy"]
        adsorbate_name = self["adsorbate_name"]
        adsite_types = self["adsite_types"]

        logging.debug(fw_spec)
        ads_structures = []
        connect_dict = {}
        reverse_connect_dict = {}
        desc_lst = []

        # going through nc atoms
        all_atomtypes = fw_spec["nc_atomic_numbers"]
        nc_structures_dict = fw_spec["nc_structures"]
        for idx, atoms_dict in enumerate(nc_structures_dict):
            logging.debug("ATOMS DICT")
            logging.debug(atoms_dict)
            atoms = ase.Atoms()
            #atoms.__dict__ = atoms_dict
            
            atoms =  atoms_dict_to_ase(atoms_dict)
            logging.debug("ASE atoms from DICT")
            logging.debug(atoms)

            # entries for adsorbate dictionary
            connect_dict[str(idx)] = {}

            # running clusgeo on cluster
            surfatoms = clusgeo.surface.get_surface_atoms(atoms)
            logging.debug(len(surfatoms))
            for adsite_type in adsite_types:
                if adsite_type == "top":
                    adsites = clusgeo.surface.get_top_sites(atoms, surfatoms)
                elif adsite_type == "bridge":
                    adsites = clusgeo.surface.get_edge_sites(atoms, surfatoms)
                elif adsite_type == "hollow":
                    adsites = clusgeo.surface.get_hollow_sites(atoms, surfatoms)
                else:
                    logging.error("adsorption site type unknown, known types are: top, bridge, hollow")
                    exit(1)
                
                # get adsorption sites for a nanocluster
                adsites_dict = adsorbate_pos_to_atoms_dict(atoms, adsites, adsorbate_name)

                # provide ids for adsorbate dictionary
                added_ids_start = len(ads_structures)
                added_ids_end = added_ids_start + len(adsites_dict)
                added_ids = list(range(added_ids_start, added_ids_end))
                connect_dict[str(idx)][adsite_type] = added_ids
                for adsorbate_id in added_ids:
                    reverse_connect_dict[str(adsorbate_id)] = str(idx)

                # save adsorbate structures
                ads_structures.extend(adsites_dict)

                # get descriptor
                desc = adsorbate_pos_to_descriptor(atoms, adsites, adsite_type, all_atomtypes, rcut = 5.0)
                for i in range(desc.shape[0]):
                    desc_lst.append(desc[i])

        descmatrix = np.array(desc_lst)
    
        # order cluster-adsorbates
        update_spec = fw_spec
        update_spec["nc_structures"] = nc_structures_dict
        update_spec["ads_structures"] = ads_structures
        update_spec["connect_dict"] = connect_dict
        update_spec["reverse_connect_dict"] = reverse_connect_dict
        update_spec["descmatrix"] = descmatrix

        # dictionary and list for filling total energies later
        update_spec["adsorbate_energies_dict"] = {}
        update_spec["adsorbate_energies_list"] = np.zeros(descmatrix.shape[0])
        
        # dictionary and list for filling structures later
        update_spec["relaxed_structure_dict"] = {}
        update_spec["relaxed_structure_list"] = list(np.zeros(descmatrix.shape[0]))

        # dictionary and list for filling reaction energies later
        update_spec["reaction_energies_list"] = np.zeros(descmatrix.shape[0])

        update_spec.pop("_category")

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
        ads_structures = fw_spec["ads_structures"]
        nc_ids = fw_spec["nc_ids"]
        descmatrix = fw_spec["descmatrix"]
        descmatrix = np.array(descmatrix)
        logging.info("DESCRIPTOR matrix attributes")
        logging.info(descmatrix.shape)
        logging.info(np.sum(descmatrix))
        fps_ranking = clusgeo.environment.rank_sites(descmatrix, K = None, idx=[], greedy = False, is_safe = True)

        update_spec = fw_spec
        update_spec["fps_ranking"] = fps_ranking
        update_spec.pop("_category")

        return FWAction(update_spec=update_spec)






def get_adsites(reference_energy = 0.0, adsorbate_name='H', adsite_types = ["top", "bridge", "hollow"],
):
    firetask1  = AdsiteCreationTask(
        reference_energy=reference_energy, 
        adsorbate_name=adsorbate_name, 
        adsite_types = adsite_types,
        )
    fw = Firework([firetask1], 
        spec={'reference_energy': reference_energy, 'adsorbate_name' : adsorbate_name, 
        'adsite_types' : adsite_types,
        '_category' : "lightweight",
}
        )
    return fw


def rank_adsites():
    firetask1  = AdsiteRankTask()
    fw = Firework([firetask1], spec = {'_category' : "medium"})
    return fw