import numpy as np
from ase.io import read

from nomadcore.atoms2nomad import ase_atoms_to_section_system

class CubeError(OSError):
    pass


def read_cube_file(backend, file_name):
    try:
        d = read(file_name, format = 'cube', full_output = True)

    except Exception as err:
        raise CubeError(err)

    data = d['data']
    atoms = d['atoms']
    origin = d['origin']
    nx, ny, nz = data.shape
    displacements = np.array([atoms.cell[i]/data.shape[i] for i in range(3)])

    system = ase_atoms_to_section_system(backend, atoms)

    singleconfig = backend.openSection('section_single_configuration_calculation')
    volumetric = backend.openSection('section_volumetric_data')

    backend.addValue('volumetric_data_nx', nx)
    backend.addValue('volumetric_data_ny', ny)
    backend.addValue('volumetric_data_nz', nz)

    backend.addArrayValues('volumetric_data_origin', origin)
    backend.addArrayValues('volumetric_data_displacements', displacements)

    backend.addValue('volumetric_data_multiplicity', 1)
    backend.addArrayValues('volumetric_data_values', data[None])
    backend.closeSection('section_volumetric_data', volumetric)

    backend.addValue('single_configuration_calculation_to_system_ref', system)
    backend.closeSection('section_single_configuration_calculation', singleconfig)
