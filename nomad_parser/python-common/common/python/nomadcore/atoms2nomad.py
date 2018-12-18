import numpy as np

from nomadcore.unit_conversion.unit_conversion import convert_unit


def ase_atoms_to_section_system(backend, atoms, new_section=True):
    """Add ASE Atoms object as metainfo to section_system.

    If new_section is True, open and close a new section_system,
    returning its gIndex."""

    if new_section:
        gIndex = backend.openSection('section_system')

    backend.addArrayValues('atom_labels',
                           np.array(atoms.get_chemical_symbols()))
    backend.addArrayValues('atom_positions',
                           convert_unit(atoms.positions, 'angstrom'))
    backend.addArrayValues('simulation_cell',
                           convert_unit(atoms.cell, 'angstrom'))
    backend.addArrayValues('configuration_periodic_dimensions',
                           np.array(atoms.pbc))

    # Return system ref if we opened it, else None:
    if new_section:
        backend.closeSection('section_system', gIndex)
        return gIndex
