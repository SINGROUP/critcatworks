from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time, pathlib, sys
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import ase, ase.io
import logging
from critcatworks.database.extdb import update_simulations_collection
from critcatworks.database.format import atoms_dict_to_ase, ase_to_atoms_dict
from critcatworks.database.extdb import fetch_simulations
import json
import numpy as np
from critcatworks.database.extdb import get_external_database, _query_id_counter_and_increment
import cluskit
import dscribe
@explicit_serialize
class NCStabilityTask(FiretaskBase):
    """ 
    Task to compare the stability of clusters using cohesive or total energy.
    Cohesive energy requires atomic total energies.

    Writes a file called cohesive_energy.json where the energies are summarized
    split by nanocluster composition.

    Args:
        atomic_energies (list) :    atomic energies provided using the same 
                                    simulation parameters

    Returns:
        FWAction : Firework action, updates fw_spec
    """
    _fw_name = 'NCStabilityTask'
    required_params = []
    optional_params = ['atomic_energies']

    def run_task(self, fw_spec):
        atomic_energies = self.get("atomic_energies", {})
        workflow_id = fw_spec.get("workflow", {"_id" : -1 }).get("_id", -1)
        logging.debug(fw_spec)
        update_spec = fw_spec
        analysis_ids = fw_spec["temp"]["analysis_ids"]
        # analysis_ids becomes calc_ids
        calc_ids = analysis_ids
        simulations = fetch_simulations(fw_spec["extdb_connect"], calc_ids)

        cohesive_energy_dct = {}
        cohesive_energy_lst = []
        for calc_id in calc_ids:
            simulation = simulations[str(calc_id)]
            print(simulation["output"])
            # get total energy
            try:
                total_energy = simulation["output"]["total_energy"]
            except:
                logging.warning("simulation" + str(calc_id) + " did not seam to have total_energy in its output.")
                continue

            atomic_numbers = simulation["atoms"]["numbers"]
            cohesive_energy = total_energy

            if cohesive_energy == None:
                logging.warning("simulation" + str(calc_id) + " has a total_energy of None.")
                continue                

            # get chemical formula
            dct = ase.utils.formula._count_symbols(atomic_numbers)
            for symbol, occurences in dct.items():
                cohesive_energy = cohesive_energy - ( atomic_energies.get(str(symbol), 0.0) * occurences )

            formula_string = ase.utils.formula.formula_metal(atomic_numbers, empirical=False)
            if formula_string not in cohesive_energy_dct:
                cohesive_energy_dct[formula_string] = {}
            cohesive_energy_dct[formula_string][str(calc_id)] = cohesive_energy
            cohesive_energy_lst.append(cohesive_energy)


        sort_ids = np.argsort(cohesive_energy_lst)
        sorted_calc_ids = np.array(calc_ids)[sort_ids]

        logging.info("cohesive_energy_dct")
        logging.info("cohesive_energy_lst")
        logging.info(cohesive_energy_dct)
        logging.info(cohesive_energy_lst)

        # cluster_stability
        with open('cohesive_energy.json', 'w') as outfile:
            json.dump(cohesive_energy_dct, outfile)

        # fireworks
        update_spec = fw_spec
        update_spec["temp"]["sorted_calc_ids"] = sorted_calc_ids
        update_spec["temp"]["cohesive_energy_dct"] = cohesive_energy_dct
        update_spec["calc_ids"] = calc_ids

        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)


@explicit_serialize
class NCGenerationTask(FiretaskBase):
    """ 
    Firetask to generate binary nanoclusters of defined size and shape.
    For each binary element combination, for each composition, n_configurations
    maximally dissimilar structures are created and uploaded to the
    simulations collection of the mongodb database.
    For further information on the generation algorithm, 
    consult the documentation of cluskit.

    The new structures are added to the simulation collection of
    the mongodb database.

    Args:
        n_initial_configurations (int) : number of initial configurations per 
                                         composition to choose from (higher number
                                         will make the grid finer)
        n_configurations (int) :    number of configurations per composition 
                                    (chosen based on maximally different 
                                    structural features)
        shape (str) : determines shape of nanoclusters. 'ico', 'octa' and 'wulff' 
        nanocluster_size (int) : determines nanocluster size. Meaning depends on shape 
        compositions (list) : each element determines the amount of atoms of type 1. 
        elements (list) : elements (str) to iterate over
        generate_pure_nanoclusters (bool) : if set to True, also pure 
                                            nanoclusters are generated
        bondlength_dct (dict) :     bond lengths to use for specific elements. 
                                    Some default bond lenghts are provided for
                                    common elements

    Returns:
        FWAction : Firework action, updates fw_spec
    """

    _fw_name = 'NCGenerationTask'
    required_params = ['n_initial_configurations', 'n_configurations',
        'shape', 'nanocluster_size', 'compositions', 'elements',]
    optional_params = ['generate_pure_nanoclusters', 'bondlength_dct']

    def run_task(self, fw_spec):
        workflow_id = fw_spec.get("workflow", {"_id" : -1 }).get("_id", -1)
        n_initial_configurations = self["n_initial_configurations"] 
        n_configurations = self["n_configurations"]
        shape = self["shape"]
        nanocluster_size = self["nanocluster_size"]
        compositions = self["compositions"]
        elements = self["elements"]
        generate_pure_nanoclusters = self["generate_pure_nanoclusters"], 
        bondlength_dct = self["bondlength_dct"]

        db = get_external_database(fw_spec["extdb_connect"])
        simulations = db['simulations']

        # generate clusters
        nanoclusters, calc_ids = self.generate(n_initial_configurations, n_configurations, 
            shape, nanocluster_size, 
            compositions, elements, 
            generate_pure_nanoclusters = generate_pure_nanoclusters, bondlength_dct = bondlength_dct,
            db = db, workflow_id = workflow_id)

        # upload all simulations at once
        simulations.insert_many(nanoclusters)

        # fireworks
        update_spec = fw_spec
        update_spec["calc_ids"] = calc_ids

        update_spec.pop("_category")
        update_spec.pop("name")
        return FWAction(update_spec=update_spec)


    def generate(self, n_initial_configurations, n_configurations, shape, 
        nanocluster_size, compositions, elements, 
        generate_pure_nanoclusters = True, bondlength_dct = {}, db = None, workflow_id = -1):
        """
        Generates binary nanoclusters. For further information, 
        consult the documentation of cluskit

        Args:
            n_initial_configurations (int) : number of initial configurations per 
                                             composition to choose from (higher number
                                             will make the grid finer)
            n_configurations (int) :    number of configurations per composition 
                                        (chosen based on maximally different 
                                        structural features)
            shape (str) : determines shape of nanoclusters. 'ico', 'octa' and 'wulff' 
            nanocluster_size (int) : determines nanocluster size. Meaning depends on shape 
            compositions (list) : each element determines the amount of atoms of type 1. 
            elements (list) : elements (str) to iterate over
            generate_pure_nanoclusters (bool) : if set to True, also pure 
                                                nanoclusters are generated
            bondlength_dct (dict) :     bond lengths to use for specific elements. 
                                        Some default bond lenghts are provided for
                                        common elements
            db (pymongo object) :       connection to the mongodb database
            workflow_id (int) :     unique identifier of the current workflow

        Returns:
            tuple : list of nanocluster ase.Atoms objects,
                    list of simulation ids
        """
        test_scaffold = cluskit.build.get_scaffold(shape, nanocluster_size, 3.0)
        NATOMS = len(test_scaffold)
        count = 0
        calc_ids = []
        nanoclusters_lst = []
        simulations = db['simulations']

        # lattice_constants fcc
        tm_dict = {'Sc': 4.6796, 'Ti': 4.1731, 'V': 4.2851, 'Cr': 4.1154, 'Mn': 1.2604, 'Fe': 4.0538, 
               'Co': 3.5456, 'Zn': 3.7687, 'Y': 5.1582, 'Zr': 4.5707, 'Nb': 4.6675, 'Mo': 4.4505, 
               'Tc': 3.8679, 'Ru': 3.8267, 'Cd': 4.2135, 'Hf': 4.5204, 'Ta': 4.6687, 'W': 4.4763, 
               'Re': 3.9046, 'Os': 3.8670, 'Hg': 4.2497}

        # manage compositions
        for i, el1 in enumerate(elements):
            for j, el2 in enumerate(elements):
                print(el1, el2)
                if j < i:
                    print("skipping double-counted combinations")
                    continue
                elif (j == i) and (generate_pure_nanoclusters == False):
                    print("skipping pure nanoclusters")        
                    continue
            

                for composition in compositions:
                    # estimate bondlength
                    bondlength_el1 = bondlength_dct.get(el1, tm_dict.get(el1, 4.0) / 1.41)
                    bondlength_el2 = bondlength_dct.get(el2, tm_dict.get(el2, 4.0) / 1.41)

                    bondlength = composition / NATOMS * bondlength_el1 + (NATOMS - composition) / NATOMS * bondlength_el2

                    # create scaffold
                    print("creating scaffold")
                    scaffold = cluskit.build.get_scaffold(shape, nanocluster_size, bondlength * 1.41)

                    atnum1 = ase.data.atomic_numbers[el1]
                    atnum2 = ase.data.atomic_numbers[el2]

                    scaffold.descriptor_setup = dscribe.descriptors.SOAP(
                        species=[atnum1, atnum2],
                        periodic=False,
                        rcut=5.0,
                        nmax=8,
                        lmax=6,
                        sparse=False,
                        average=True
                    )
                    # run cluster generator
                    print("cluster generator")
                    if i == j:
                        cluster_lst = scaffold.get_unique_clusters_in_range(typeA = atnum1, typeB = atnum2, 
                            ntypeB = NATOMS - composition, n_clus = 1)
                    else:
                        cluster_lst = scaffold.get_unique_clusters_in_range(typeA = atnum1, typeB = atnum2, 
                            ntypeB = NATOMS - composition, n_clus = n_initial_configurations)
                        cluster_lst = cluster_lst[:n_configurations]

                    count += n_configurations
                    for count, cluster in enumerate(cluster_lst):
                        name = el1 + str(composition) + el2 + str(NATOMS - composition) + "_" + str(count) + ".xyz"
                        # simulation format
                        simulation_id = _query_id_counter_and_increment('simulations', db)
                        atoms = cluster.ase_object
                        nanocluster_atom_ids = list(range(len(atoms)))
                        atoms_dict = ase_to_atoms_dict(atoms)

                        nanocluster = {"atom_ids" : nanocluster_atom_ids, "reference_id" : simulation_id}

                        dct = {"_id" : simulation_id, "atoms" : atoms_dict, 
                            "source_id" : -1, "workflow_id" : workflow_id, 
                            "nanoclusters" : [nanocluster], "adsorbates" : [], "substrates" : [], 
                            "operations" : [], "inp" : cluster.info,
                            }

                        nanoclusters_lst.append(dct)
                        calc_ids.append(simulation_id)

                    if i == j:
                        break
        return nanoclusters_lst, calc_ids


def compare_nanoclusters(atomic_energies = {}):    
    """ 
    Firework to compare the stability of clusters using cohesive or total energy.
    Cohesive energy requires atomic total energies.

    Writes a file called cohesive_energy.json where the energies are summarized
    split by nanocluster composition.

    Args:
        atomic_energies (list) :    atomic energies provided using the same 
                                    simulation parameters

    Returns:
        Firework : Firework NCStabilityWork
    """
    firetask1  = NCStabilityTask(atomic_energies = atomic_energies)
    dct = {'_category' : "lightweight", 'name' : 'NCStabilityTask'}
    fw = Firework([firetask1], spec=dct,
             name = 'NCStabilityWork')
    return fw


def generate_nanoclusters(n_initial_configurations, n_configurations, 
    shape, nanocluster_size, compositions, elements, 
    generate_pure_nanoclusters = True, bondlength_dct = {}):
    """
    Firework to generate binary nanoclusters of defined size and shape.
    For each binary element combination, for each composition, n_configurations
    maximally dissimilar structures are created and uploaded to the
    simulations collection of the mongodb database.
    For further information on the generation algorithm, 
    consult the documentation of cluskit.

    The new structures are added to the simulation collection of
    the mongodb database.

    Args:
        n_initial_configurations (int) : number of initial configurations per 
                                         composition to choose from (higher number
                                         will make the grid finer)
        n_configurations (int) :    number of configurations per composition 
                                    (chosen based on maximally different 
                                    structural features)
        shape (str) : determines shape of nanoclusters. 'ico', 'octa' and 'wulff' 
        nanocluster_size (int) : determines nanocluster size. Meaning depends on shape 
        compositions (list) : each element determines the amount of atoms of type 1. 
        elements (list) : elements (str) to iterate over
        generate_pure_nanoclusters (bool) : if set to True, also pure 
                                            nanoclusters are generated
        bondlength_dct (dict) :     bond lengths to use for specific elements. 
                                    Some default bond lenghts are provided for
                                    common elements

    Returns:
        Firework : Firework NCGenerationWork
    """
    firetask1  = NCGenerationTask(n_initial_configurations = n_initial_configurations, 
        n_configurations = n_configurations, shape = shape, nanocluster_size = nanocluster_size,
        compositions = compositions, elements = elements,
        generate_pure_nanoclusters = True, bondlength_dct = {})
    dct = {'_category' : "large", 'name' : 'NCGenerationTask'}
    fw = Firework([firetask1], spec=dct,
             name = 'NCGenerationWork')
    return fw
