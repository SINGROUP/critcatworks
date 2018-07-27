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

def atoms_dict_to_ase(atoms_dict):
    calculator = atoms_dict['_calc']
    cell = atoms_dict['_cell']
    celldisp = atoms_dict['_celldisp']
    constraint = atoms_dict['_constraints']
    pbc = atoms_dict['_pbc']
    numbers = atoms_dict['arrays']['numbers']
    positions = atoms_dict['arrays']['positions']
    info = atoms_dict['info']

    atoms = ase.Atoms(numbers=numbers,
        positions=positions,
        calculator=calculator,
        cell = cell,
        celldisp = celldisp,
        constraint = constraint,
        pbc = pbc,
        info = info)

    return atoms

"""
    {'_calc': None,
 '_cell': array([[0., 0., 0.],
       [0., 0., 0.],
       [0., 0., 0.]]),
 '_celldisp': array([[0.],
       [0.],
       [0.]]),
 '_constraints': [],
 '_pbc': array([False, False, False]),
 'arrays': {'numbers': array([8, 1, 1]),
            'positions': array([[ 0.      ,  0.      ,  0.119262],
       [ 0.      ,  0.763239, -0.477047],
       [ 0.      , -0.763239, -0.477047]])},
 'info': {}}
"""

@explicit_serialize
class AdsiteCreationTask(FiretaskBase):
    """ 
    Task to determine adsorption site structures.

    Args:
        adsorbate_energy (float) :  Total energy of the free adsorbate 
                                    as a reference energy to 
                                    determine the adsorption energy 
                                    as in
                                    dE = dE(all) - dE(nc) - dE(adsorbate)
        adsorbate_name  (str) : Adsorbate atom name to be placed
                                on all sites found.
    """

    _fw_name = 'AdsiteCreationTask'
    required_params = ['adsorbate_energy', 'adsorbate_name']
    optional_params = []

    def run_task(self, fw_spec):
        adsorbate_energy = self["adsorbate_energy"]
        adsorbate_name = self["adsorbate_name"]

        pp(fw_spec)

        # going through nc atoms
        nc_structures_dict = fw_spec["nc_structures"]
        for atoms_dict in nc_structures_dict:
            print("ATOMS DICT")
            print(atoms_dict)
            atoms = ase.Atoms()
            #atoms.__dict__ = atoms_dict
            
            atoms =  atoms_dict_to_ase(atoms_dict)
            print("ASE atoms from DICT")
            print(atoms)
            

            # running clusgeo on cluster
            surfatoms = clusgeo.surface.get_surface_atoms(atoms)
            top_adsites = clusgeo.surface.get_top_sites(atoms, surfatoms)
            edge_adsites = clusgeo.surface.get_edge_sites(atoms, surfatoms)
            hollow_adsites = clusgeo.surface.get_hollow_sites(atoms, surfatoms)

        print(len(surfatoms))
        print(len(top_adsites))
        print(len(edge_adsites))
        print(len(hollow_adsites))

        update_spec = fw_spec
        update_spec["nc_structures"] = nc_structures_dict

        return FWAction(update_spec=update_spec)




def get_adsites(adsorbate_energy = 0.0, adsorbate_name='H'):
    firetask1  = AdsiteCreationTask(adsorbate_energy=adsorbate_energy, 
        adsorbate_name=adsorbate_name, 
        )
    fw = Firework([firetask1], 
        spec={'adsorbate_energy': adsorbate_energy, 'adsorbate_name' : adsorbate_name}
        )
    return fw
