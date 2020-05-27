# Copyright 2015-2018 Lauri Himanen, Fawzi Mohamed, Ankit Kariryaa
# 
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

from __future__ import absolute_import
from builtins import str
import re
import numpy as np
import logging
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.caching_backend import CachingLevel
from nomadcore.unit_conversion.unit_conversion import convert_unit
from nomadcore.baseclasses import CommonParser
from .inputparser import CP2KInputParser
from collections import defaultdict
logger = logging.getLogger("nomad")


class CP2KCommonParser(CommonParser):
    """
    This class is used to store and instantiate common parts of the
    hierarchical SimpleMatcher structure used in the parsing of a CP2K
    output file.
    """
    def __init__(self, parser_context):
        super(CP2KCommonParser, self).__init__(parser_context)
        self.section_method_index = None
        self.section_system_index = None
        self.test_electronic_structure_method = "DFT"
        self.basis_to_kind_mapping = []
        self.atom_kind_info = defaultdict(dict)  # Map from kind number to kind information
        self.basis_set_info = defaultdict(dict)  # Map from kind number to basis set information

        #=======================================================================
        # Cache levels
        self.caching_levels = {
            'x_cp2k_atoms': CachingLevel.ForwardAndCache,
            'section_XC_functionals': CachingLevel.ForwardAndCache,
            'self_interaction_correction_method': CachingLevel.Cache,
            'x_cp2k_section_program_information': CachingLevel.ForwardAndCache,
            'x_cp2k_section_quickstep_settings': CachingLevel.ForwardAndCache,
            'x_cp2k_section_atomic_kind': CachingLevel.ForwardAndCache,
            'x_cp2k_section_kind_basis_set': CachingLevel.ForwardAndCache,
        }

        #=======================================================================
        # Globally cached values
        self.cache_service.add("simulation_cell", single=False, update=False)
        self.cache_service.add("lattice_vectors", single=False, update=False)
        self.cache_service.add("number_of_scf_iterations", 0)
        self.cache_service.add("atom_positions", single=False, update=True)
        self.cache_service.add("atom_labels", single=False, update=False)
        self.cache_service.add("number_of_atoms", single=False, update=False)
        self.cache_service.add("basis_set_kind", single=False, update=False)
        self.cache_service.add("basis_set_name", single=False, update=False)
        self.cache_service.add("basis_set_planewave_cutoff", update=False)
        self.cache_service.add("mapping_section_basis_set_cell_dependent", single=False, update=False)
        self.cache_service.add("map_kind_to_basis", single=False, update=False)
        self.cache_service.add("map_index_to_kind", single=False, update=False)
        self.cache_service.add("map_kind_number_to_basis_ref", single=False, update=False)
        self.cache_service.add("electron_density_filename", single=False, update=True)

    #===========================================================================
    # SimpleMatchers

    # SimpleMatcher for the header that is common to all run types
    def header(self):
        return SM( " DBCSR\| Multiplication driver",
            forwardMatch=True,
            subMatchers=[
                SM( " DBCSR\| Multiplication driver",
                    forwardMatch=True,
                    sections=['x_cp2k_section_dbcsr'],
                    subMatchers=[
                        SM( " DBCSR\| Multiplication driver\s+(?P<x_cp2k_dbcsr_multiplication_driver>{})".format(self.regexs.word)),
                        SM( " DBCSR\| Multrec recursion limit\s+(?P<x_cp2k_dbcsr_multrec_recursion_limit>{})".format(self.regexs.int)),
                        SM( " DBCSR\| Multiplication stack size\s+(?P<x_cp2k_dbcsr_multiplication_stack_size>{})".format(self.regexs.int)),
                        SM( " DBCSR\| Multiplication size stacks\s+(?P<x_cp2k_dbcsr_multiplication_size_stacks>{})".format(self.regexs.int)),
                        SM( " DBCSR\| Use subcommunicators\s+(?P<x_cp2k_dbcsr_use_subcommunicators>{})".format(self.regexs.letter)),
                        SM( " DBCSR\| Use MPI combined types\s+(?P<x_cp2k_dbcsr_use_mpi_combined_types>{})".format(self.regexs.letter)),
                        SM( " DBCSR\| Use MPI memory allocation\s+(?P<x_cp2k_dbcsr_use_mpi_memory_allocation>{})".format(self.regexs.letter)),
                        SM( " DBCSR\| Use Communication thread\s+(?P<x_cp2k_dbcsr_use_communication_thread>{})".format(self.regexs.letter)),
                        SM( " DBCSR\| Communication thread load\s+(?P<x_cp2k_dbcsr_communication_thread_load>{})".format(self.regexs.int)),
                    ]
                ),
                SM( "  **** **** ******  **  PROGRAM STARTED AT".replace("*", "\*"),
                    forwardMatch=True,
                    sections=['x_cp2k_section_startinformation'],
                    subMatchers=[
                        SM( "  **** **** ******  **  PROGRAM STARTED AT\s+(?P<x_cp2k_start_time>{})".replace("*", "\*").format(self.regexs.eol)),
                        SM( " ***** ** ***  *** **   PROGRAM STARTED ON\s+(?P<x_cp2k_start_host>{})".replace("*", "\*").format(self.regexs.word)),
                        SM( " **    ****   ******    PROGRAM STARTED BY\s+(?P<x_cp2k_start_user>{})".replace("*", "\*").format(self.regexs.word)),
                        SM( " ***** **    ** ** **   PROGRAM PROCESS ID\s+(?P<x_cp2k_start_id>{})".replace("*", "\*").format(self.regexs.int)),
                        SM( "  **** **  *******  **  PROGRAM STARTED IN".replace("*", "\*"),
                            forwardMatch=True,
                            adHoc=self.adHoc_run_dir("x_cp2k_start_path"),
                        )
                    ]
                ),
                SM( " CP2K\| version string:",
                    sections=['x_cp2k_section_program_information'],
                    forwardMatch=True,
                    subMatchers=[
                        SM( " CP2K\| version string:\s+(?P<program_version>{})".format(self.regexs.eol)),
                        SM( " CP2K\| source code revision number:\s+svn:(?P<x_cp2k_svn_revision>\d+)"),
                        SM( " CP2K\| is freely available from{}".format(self.regexs.eol)),
                        SM( " CP2K\| Program compiled at\s+(?P<x_cp2k_program_compilation_datetime>{})".format(self.regexs.eol)),
                        SM( " CP2K\| Program compiled on\s+(?P<program_compilation_host>{})".format(self.regexs.eol)),
                        SM( " CP2K\| Program compiled for{}".format(self.regexs.eol)),
                        SM( " CP2K\| Input file name\s+(?P<x_cp2k_input_filename>{})".format(self.regexs.eol)),
                    ]
                ),
                SM( " GLOBAL\|",
                    sections=['x_cp2k_section_global_settings'],
                    subMatchers=[
                        SM( " GLOBAL\| Force Environment number"),
                        SM( " GLOBAL\| Basis set file name\s+(?P<x_cp2k_basis_set_filename>{})".format(self.regexs.eol)),
                        SM( " GLOBAL\| Geminal file name\s+(?P<x_cp2k_geminal_filename>{})".format(self.regexs.eol)),
                        SM( " GLOBAL\| Potential file name\s+(?P<x_cp2k_potential_filename>{})".format(self.regexs.eol)),
                        SM( " GLOBAL\| MM Potential file name\s+(?P<x_cp2k_mm_potential_filename>{})".format(self.regexs.eol)),
                        SM( " GLOBAL\| Coordinate file name\s+(?P<x_cp2k_coordinate_filename>{})".format(self.regexs.eol)),
                        SM( " GLOBAL\| Method name\s+(?P<x_cp2k_method_name>{})".format(self.regexs.eol)),
                        SM( " GLOBAL\| Project name"),
                        SM( " GLOBAL\| Preferred FFT library\s+(?P<x_cp2k_preferred_fft_library>{})".format(self.regexs.eol)),
                        SM( " GLOBAL\| Preferred diagonalization lib.\s+(?P<x_cp2k_preferred_diagonalization_library>{})".format(self.regexs.eol)),
                        SM( " GLOBAL\| Run type\s+(?P<x_cp2k_run_type>{})".format(self.regexs.eol)),
                        SM( " GLOBAL\| All-to-all communication in single precision"),
                        SM( " GLOBAL\| FFTs using library dependent lengths"),
                        SM( " GLOBAL\| Global print level"),
                        SM( " GLOBAL\| Total number of message passing processes"),
                        SM( " GLOBAL\| Number of threads for this process"),
                        SM( " GLOBAL\| This output is from process"),
                    ],
                ),
                SM( " CELL\|",
                    adHoc=self.adHoc_x_cp2k_section_cell,
                ),
            ]
        )

    # SimpleMatcher for the restart information
    def restart(self):
        return SM( re.escape(" *                            RESTART INFORMATION                              *"),
            sections=["x_cp2k_section_restart_information"],
            subMatchers=[
                SM( re.escape(" *******************************************************************************")),
                SM( re.escape(" *                                                                             *")),
                SM( re.escape(" *    RESTART FILE NAME: (?P<x_cp2k_restart_file_name>{})\s+*".format(self.regexs.word))),
                SM( re.escape(" *                                                                             *")),
                SM( re.escape(" * RESTARTED QUANTITIES:                                                       *")),
                SM( re.escape(" *                       (?P<x_cp2k_restarted_quantity>{})\s+*".format(self.regexs.word)),
                    repeats=True,
                ),
                SM( re.escape(" *******************************************************************************"))
            ]
        )

    # SimpleMatcher for the footer that is common to all run types
    def footer(self):
        return SM( " -                                DBCSR STATISTICS                             -",
            forwardMatch=True,
            subMatchers=[
                SM( re.escape("  **** **** ******  **  PROGRAM ENDED AT"),
                    forwardMatch=True,
                    sections=['x_cp2k_section_end_information'],
                    subMatchers=[
                        SM( "  **** **** ******  **  PROGRAM ENDED AT\s+(?P<x_cp2k_end_time>{})".replace("*", "\*").format(self.regexs.eol)),
                        SM( " ***** ** ***  *** **   PROGRAM RAN ON\s+(?P<x_cp2k_end_host>{})".replace("*", "\*").format(self.regexs.word)),
                        SM( " **    ****   ******    PROGRAM RAN BY\s+(?P<x_cp2k_end_user>{})".replace("*", "\*").format(self.regexs.word)),
                        SM( " ***** **    ** ** **   PROGRAM PROCESS ID\s+(?P<x_cp2k_end_id>{})".replace("*", "\*").format(self.regexs.int)),
                        SM( "  **** **  *******  **  PROGRAM STOPPED IN".replace("*", "\*"),
                            forwardMatch=True,
                            adHoc=self.adHoc_run_dir("x_cp2k_end_path"),
                        )
                    ]
                ),
            ]
        )

    # SimpleMatcher for an SCF wavefunction optimization
    def quickstep_calculation(self):
        return SM( " SCF WAVEFUNCTION OPTIMIZATION",
            sections=["x_cp2k_section_quickstep_calculation"],
            subMatchers=[
                SM( r"  Trace\(PS\):",
                    sections=["x_cp2k_section_scf_iteration"],
                    repeats=True,
                    subMatchers=[
                        SM( r"  Exchange-correlation energy:\s+(?P<x_cp2k_energy_XC_scf_iteration__hartree>{})".format(self.regexs.float)),
                        SM( r"\s+\d+\s+\S+\s+{0}\s+{0}\s+{0}\s+(?P<x_cp2k_energy_total_scf_iteration__hartree>{0})\s+(?P<x_cp2k_energy_change_scf_iteration__hartree>{0})".format(self.regexs.float)),
                    ]
                ),
                SM( r"  \*\*\* SCF run converged in\s+(\d+) steps \*\*\*",
                    adHoc=self.adHoc_single_point_converged
                ),
                SM( r"  \*\*\* SCF run NOT converged \*\*\*",
                    adHoc=self.adHoc_single_point_not_converged
                ),
                SM( r" The electron density is written in cube file format to the file:",
                    subMatchers=[
                        SM(""),
                        SM(
                            "\s+(.+\.cube)",
                            startReAction=self.startReAction_save_cube_filename,
                        ),
                    ]
                ),
                SM( r"  Electronic kinetic energy:\s+(?P<x_cp2k_electronic_kinetic_energy__hartree>{})".format(self.regexs.float)),
                SM( r" **************************** NUMERICAL STRESS ********************************".replace("*", "\*"),
                    # endReStr=" **************************** NUMERICAL STRESS END *****************************".replace("*", "\*"),
                    adHoc=self.adHoc_stress_calculation,
                ),
                SM( r" ENERGY\| Total FORCE_EVAL \( \w+ \) energy \(a\.u\.\):\s+(?P<x_cp2k_energy_total__hartree>{0})".format(self.regexs.float),
                ),
                SM( r" ATOMIC FORCES in \[a\.u\.\]"),
                SM( r" # Atom   Kind   Element          X              Y              Z",
                    adHoc=self.adHoc_atom_forces,
                ),
                SM( r" (?:NUMERICAL )?STRESS TENSOR \[GPa\]",
                    sections=["x_cp2k_section_stress_tensor"],
                    subMatchers=[
                        SM( r"\s+X\s+Y\s+Z",
                            adHoc=self.adHoc_stress_tensor,
                        ),
                        SM( "  1/3 Trace\(stress tensor\):\s+(?P<x_cp2k_stress_tensor_one_third_of_trace__GPa>{})".format(self.regexs.float)),
                        SM( "  Det\(stress tensor\)\s+:\s+(?P<x_cp2k_stress_tensor_determinant__GPa3>{})".format(self.regexs.float)),
                        SM( " EIGENVECTORS AND EIGENVALUES OF THE STRESS TENSOR",
                            adHoc=self.adHoc_stress_tensor_eigenpairs),
                    ]
                )
            ]
        )

    # SimpleMatcher the stuff that is done to initialize a quickstep calculation
    def quickstep_header(self):
        return SM( re.escape(" **                                                ... make the atoms dance   **"),
            forwardMatch=True,
            sections=["x_cp2k_section_quickstep_settings"],
            subMatchers=[
                SM( " DFT\|",
                    forwardMatch=True,
                    subMatchers=[
                        SM( " DFT\| Spin restricted Kohn-Sham (RKS) calculation\s+(?P<x_cp2k_spin_restriction>{})".format(self.regexs.word)),
                        SM( " DFT\| Multiplicity\s+(?P<spin_target_multiplicity>{})".format(self.regexs.int)),
                        SM( " DFT\| Number of spin states\s+(?P<number_of_spin_channels>{})".format(self.regexs.int)),
                        SM( " DFT\| Charge\s+(?P<total_charge>{})".format(self.regexs.int)),
                        SM( " DFT\| Self-interaction correction \(SIC\)\s+(?P<self_interaction_correction_method>[^\n]+)"),
                    ],
                ),
                SM( " vdW POTENTIAL\|\s+",
                    forwardMatch=True,
                    sections=["x_cp2k_section_vdw_settings"],
                    subMatchers=[
                        SM( " vdW POTENTIAL\|\s+(?P<x_cp2k_vdw_type>{})".format(self.regexs.eol)),
                        SM( " vdW POTENTIAL\|\s+Potential Form:\s+(?P<x_cp2k_vdw_name>{})".format(self.regexs.eol)),
                        SM( " vdW POTENTIAL\|\s+BJ Damping:\s+(?P<x_cp2k_vdw_bj_damping_name>{})".format(self.regexs.eol)),
                        SM( " vdW POTENTIAL\|\s+Cutoff Radius \[Bohr\]:\s+(?P<x_cp2k_vdw_cutoff_radius>{})".format(self.regexs.float)),
                        SM( " vdW POTENTIAL\|\+sScaling Factor:s+(?P<x_cp2k_vdw_scaling_factor>{})".format(self.regexs.float),
                            sections=["x_cp2k_section_vdw_d2_settings"],
                            subMatchers=[
                                SM( " vdW POTENTIAL\|\s+Exp Prefactor for Damping:\s+(?P<x_cp2k_vdw_damping_factor>{})".format(self.regexs.float)),
                                SM( " vdW PARAMETER\|\s+Atom=(?P<x_cp2k_vdw_parameter_element_name>{})\s+C6\[J*nm^6*mol^-1\]=\s+(?P<x_cp2k_vdw_parameter_c6>{})\s+r\(vdW\)\[A\]=\s+(?P<x_cp2k_vdw_parameter_radius>{})".format(self.regexs.word, self.regexs.float, self.regexs.float),
                                    sections=["x_cp2k_section_vdw_element_settings"],
                                    repeats=True,
                                ),
                            ],
                        ),
                        SM( " vdW POTENTIAL\|\s+s6 Scaling Factor:\s+(?P<x_cp2k_vdw_s6_scaling_factor>{})".format(self.regexs.float),
                            sections=["x_cp2k_section_vdw_d3_settings"],
                            subMatchers=[
                                SM( " vdW POTENTIAL\|\s+sr6 Scaling Factor:\s+(?P<x_cp2k_vdw_sr6_scaling_factor>{})".format(self.regexs.float)),
                                SM( " vdW POTENTIAL\|\s+s8 Scaling Factor:\s+(?P<x_cp2k_vdw_s8_scaling_factor>{})".format(self.regexs.float)),
                                SM( " vdW POTENTIAL\|\s+Cutoff for CN calculation:\s+(?P<x_cp2k_vdw_cn_cutoff>{})".format(self.regexs.float)),
                            ],
                        )
                    ]
                ),
                SM( " DFT\+U\|",
                    adHoc=self.adHoc_dft_plus_u,
                ),
                SM( " QS\|",
                    forwardMatch=True,
                    subMatchers=[
                        SM( " QS\| Method:\s+(?P<x_cp2k_quickstep_method>{})".format(self.regexs.word)),
                        SM( " QS\| Density plane wave grid type\s+{}".format(self.regexs.eol)),
                        SM( " QS\| Number of grid levels:\s+{}".format(self.regexs.int)),
                        SM( " QS\| Density cutoff \[a\.u\.\]:\s+(?P<x_cp2k_planewave_cutoff>{})".format(self.regexs.float)),
                        SM( " QS\| Multi grid cutoff \[a\.u\.\]: 1\) grid level\s+{}".format(self.regexs.float)),
                        SM( " QS\|                           2\) grid level\s+{}".format(self.regexs.float)),
                        SM( " QS\|                           3\) grid level\s+{}".format(self.regexs.float)),
                        SM( " QS\|                           4\) grid level\s+{}".format(self.regexs.float)),
                        SM( " QS\| Grid level progression factor:\s+{}".format(self.regexs.float)),
                        SM( " QS\| Relative density cutoff \[a\.u\.\]:".format(self.regexs.float)),
                        SM( " QS\| Consistent realspace mapping and integration"),
                        SM( " QS\| Interaction thresholds: eps_pgf_orb:\s+{}".format(self.regexs.float)),
                        SM( " QS\|                         eps_filter_matrix:\s+{}".format(self.regexs.float)),
                        SM( " QS\|                         eps_core_charge:\s+{}".format(self.regexs.float)),
                        SM( " QS\|                         eps_rho_gspace:\s+{}".format(self.regexs.float)),
                        SM( " QS\|                         eps_rho_rspace:\s+{}".format(self.regexs.float)),
                        SM( " QS\|                         eps_gvg_rspace:\s+{}".format(self.regexs.float)),
                        SM( " QS\|                         eps_ppl:\s+{}".format(self.regexs.float)),
                        SM( " QS\|                         eps_ppnl:\s+{}".format(self.regexs.float)),
                    ],
                ),
                SM( " ATOMIC KIND INFORMATION",
                    sections=["x_cp2k_section_atomic_kinds"],
                    subMatchers=[
                        SM( "\s+(?P<x_cp2k_kind_number>{0})\. Atomic kind: (?P<x_cp2k_kind_label>{1})\s+Number of atoms:\s+(?P<x_cp2k_kind_number_of_atoms>{1})".format(self.regexs.int, self.regexs.word),
                            repeats=True,
                            sections=["x_cp2k_section_atomic_kind", "x_cp2k_section_kind_basis_set"],
                            subMatchers=[
                                SM( "     Orbital Basis Set\s+(?P<x_cp2k_kind_basis_set_name>{})".format(self.regexs.word)),
                                SM( "       Number of orbital shell sets:\s+(?P<x_cp2k_basis_set_number_of_orbital_shell_sets>{})".format(self.regexs.int)),
                                SM( "       Number of orbital shells:\s+(?P<x_cp2k_basis_set_number_of_orbital_shells>{})".format(self.regexs.int)),
                                SM( "       Number of primitive Cartesian functions:\s+(?P<x_cp2k_basis_set_number_of_primitive_cartesian_functions>{})".format(self.regexs.int)),
                                SM( "       Number of Cartesian basis functions:\s+(?P<x_cp2k_basis_set_number_of_cartesian_basis_functions>{})".format(self.regexs.int)),
                                SM( "       Number of spherical basis functions:\s+(?P<x_cp2k_basis_set_number_of_spherical_basis_functions>{})".format(self.regexs.int)),
                                SM( "       Norm type:\s+(?P<x_cp2k_basis_set_norm_type>{})".format(self.regexs.int)),
                            ]
                        )
                    ]
                ),
                SM( "  Total number of",
                    forwardMatch=True,
                    sections=["x_cp2k_section_total_numbers"],
                    subMatchers=[
                        SM( "  Total number of            - Atomic kinds:\s+(?P<x_cp2k_atomic_kinds>\d+)"),
                        SM( "\s+- Atoms:\s+(?P<x_cp2k_atoms>\d+)"),
                        SM( "\s+- Shell sets:\s+(?P<x_cp2k_shell_sets>\d+)"),
                        SM( "\s+- Shells:\s+(?P<x_cp2k_shells>\d+)"),
                        SM( "\s+- Primitive Cartesian functions:\s+(?P<x_cp2k_primitive_cartesian_functions>\d+)"),
                        SM( "\s+- Cartesian basis functions:\s+(?P<x_cp2k_cartesian_basis_functions>\d+)"),
                        SM( "\s+- Spherical basis functions:\s+(?P<x_cp2k_spherical_basis_functions>\d+)"),
                    ]
                ),
                SM( " Maximum angular momentum of",
                    forwardMatch=True,
                    sections=["x_cp2k_section_maximum_angular_momentum"],
                    subMatchers=[
                        SM( "  Maximum angular momentum of- Orbital basis functions::\s+(?P<x_cp2k_orbital_basis_functions>\d+)"),
                        SM( "\s+- Local part of the GTH pseudopotential:\s+(?P<x_cp2k_local_part_of_gth_pseudopotential>\d+)"),
                        SM( "\s+- Non-local part of the GTH pseudopotential:\s+(?P<x_cp2k_non_local_part_of_gth_pseudopotential>\d+)"),
                    ]
                ),
                SM( " MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom",
                    forwardMatch=True,
                    subMatchers=[
                        SM( " MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom",
                            adHoc=self.adHoc_x_cp2k_section_quickstep_atom_information(),
                        )
                    ]
                ),
                SM( " SCF PARAMETERS",
                    forwardMatch=True,
                    subMatchers=[
                        SM( " SCF PARAMETERS         Density guess:\s+{}".format(self.regexs.eol)),
                        SM( "                        max_scf:\s+(?P<scf_max_iteration>{})".format(self.regexs.int)),
                        SM( "                        max_scf_history:\s+{}".format(self.regexs.int)),
                        SM( "                        max_diis:\s+{}".format(self.regexs.int)),
                        SM( "                        eps_scf:\s+(?P<scf_threshold_energy_change__hartree>{})".format(self.regexs.float)),
                    ]
                ),
                SM( " MP2\|",
                    adHoc=self.adHoc_mp2
                ),
                SM( " RI-RPA\|",
                    adHoc=self.adHoc_rpa
                ),
            ]
        )

    #===========================================================================
    # onClose triggers
    def onClose_x_cp2k_section_total_numbers(self, backend, gIndex, section):
        """Keep track of how many SCF iteration are made."""
        number_of_atoms = section.get_latest_value("x_cp2k_atoms")
        if number_of_atoms is not None:
            self.cache_service["number_of_atoms"] = number_of_atoms

    # def onClose_x_cp2k_section_quickstep_calculation(self, backend, gIndex, section):
        # print "quickstep CLOSED"

    # def onClose_x_cp2k_section_geometry_optimization_step(self, backend, gIndex, section):
        # print "Optimisation step CLOSED"

    def onClose_section_method(self, backend, gIndex, section):
        """When all the functional definitions have been gathered, matches them
        with the nomad correspondents and combines into one single string which
        is put into the backend.
        """
        self.section_method_index = gIndex

        # Transform the CP2K self-interaction correction string to the NOMAD
        # correspondent, and push directly to the superBackend to avoid caching
        try:
            sic_cp2k = section.get_latest_value("self_interaction_correction_method")
            if sic_cp2k is not None:
                sic_map = {
                    "NO": "",
                    "AD SIC": "SIC_AD",
                    "Explicit Orbital SIC": "SIC_EXPLICIT_ORBITALS",
                    "SPZ/MAURI SIC": "SIC_MAURI_SPZ",
                    "US/MAURI SIC": "SIC_MAURI_US",
                }
                sic_nomad = sic_map.get(sic_cp2k)
                if sic_nomad is not None:
                    backend.superBackend.addValue('self_interaction_correction_method', sic_nomad)
                else:
                    logger.warning("Unknown self-interaction correction method used.")
        except:
            pass

    def onClose_section_run(self, backend, gIndex, section):
        backend.addValue("program_name", "CP2K")

    def onClose_x_cp2k_section_quickstep_settings(self, backend, gIndex, section):
        backend.addValue("program_basis_set_type", "gaussians")
        backend.addValue("electronic_structure_method", self.test_electronic_structure_method)

        # Collect the atomic kind information and push it to backend
        kind_ids = {}
        for kind_number, info in self.atom_kind_info.items():
            kindID = backend.openSection("section_method_atom_kind")
            kind_ids[kind_number] = kindID
            label = info["label"]
            atom_number = info.get("element_number")
            if atom_number is not None:
                backend.addValue("method_atom_kind_atom_number", atom_number)
            backend.addValue("method_atom_kind_label", label)
            backend.closeSection("section_method_atom_kind", kindID)

        # Cell dependent basis information
        ryd_cutoff = section.get_latest_value("x_cp2k_planewave_cutoff")
        if ryd_cutoff is not None:
            gid = backend.openSection("section_basis_set_cell_dependent")
            self.cache_service["mapping_section_basis_set_cell_dependent"] = gid
            ha_cutoff = convert_unit(2*ryd_cutoff, "hartree")
            backend.addValue("basis_set_planewave_cutoff", ha_cutoff)
            self.cache_service["basis_set_planewave_cutoff"] = ryd_cutoff
            backend.closeSection("section_basis_set_cell_dependent", gid)

        # Atom centered basis set information
        basis_ids = {}
        map_kind_number_to_basis_ref = {}
        for kind_number, info in self.basis_set_info.items():
            basis_section_id = backend.openSection("section_basis_set_atom_centered")
            basis_ids[kind_number] = basis_section_id
            map_kind_number_to_basis_ref[kind_number] = basis_section_id
            name = info["name"]
            atom_number = info.get("element_number")
            if atom_number is not None:
                backend.addValue("basis_set_atom_number", atom_number)
            backend.addValue("basis_set_atom_centered_short_name", name)
            backend.closeSection("section_basis_set_atom_centered", basis_section_id)
        self.cache_service["map_kind_number_to_basis_ref"] = map_kind_number_to_basis_ref

        # Add the basis infomation to section_method
        mapping = []
        dict_map = {}
        for kind_number, basis_id in basis_ids.items():
            kind_id = kind_ids[kind_number]
            mapping.append((basis_id, kind_id))
            dict_map[kind_id] = basis_id
        method_basis_id = backend.openSection("section_method_basis_set")
        if mapping:
            mapping = np.array(mapping)
            self.cache_service["map_kind_to_basis"] = dict_map
            backend.addArrayValues("mapping_section_method_basis_set_atom_centered", np.array(mapping))
        backend.addValue("method_basis_set_kind", "wavefunction")
        self.cache_service.addValue("mapping_section_method_basis_set_cell_associated")
        backend.addValue("number_of_basis_sets_atom_centered", len(self.basis_set_info))
        backend.closeSection("section_method_basis_set", method_basis_id)

    def onClose_x_cp2k_section_vdw_settings(self, backend, gIndex, section):
        """Figures out the common metainfo for vdw from the CP2K specific
        settings.
        """
        vdw_name = section["x_cp2k_vdw_name"][0]
        name_map = {
            "S. Grimme, JCC 27: 1787 (2006)": "G06",
            "S. Grimme et al, JCP 132: 154104 (2010)": "G10",
        }
        nomad_vdw_name = name_map.get(vdw_name)
        if nomad_vdw_name is not None:
            backend.addValue("van_der_Waals_method", nomad_vdw_name)

    def onClose_x_cp2k_section_atomic_kind(self, backend, gIndex, section):
        # basisID = backend.openSection("section_basis_set_atom_centered")

        # Save the kind labels. These wil be connected to atomic numbers later
        # on when the atomic numbers are listed in the atomic positions.
        kind_number = int(section.get_latest_value("x_cp2k_kind_number"))
        kind_label = section.get_latest_value("x_cp2k_kind_label")
        self.atom_kind_info[kind_number]["label"] = kind_label

        # Save all the basis set information for later use. They will be pushed
        # later when an atom number can be associated with the basis.
        basis_set_name = section.get_latest_value(["x_cp2k_section_kind_basis_set", "x_cp2k_kind_basis_set_name"])
        basis_info = self.basis_set_info[kind_number]
        basis_info["name"] = basis_set_name

    def onClose_x_cp2k_section_atomic_kinds(self, backend, gIndex, section):
        # Store the name and kind of the basis set for later use (stored inside
        # single_configuration_calculation).
        atomic_kinds = section["x_cp2k_section_atomic_kind"]
        long_basis_name = []
        for kind in atomic_kinds:
            kind_basis = kind["x_cp2k_section_kind_basis_set"][0]
            basis_name = kind_basis["x_cp2k_kind_basis_set_name"][0]
            long_basis_name.append(basis_name)
        self.cache_service["basis_set_kind"] = "wavefunction"
        self.cache_service["basis_set_name"] = "_".join(long_basis_name)

    def onClose_x_cp2k_section_program_information(self, backend, gIndex, section):
        input_file = section.get_latest_value("x_cp2k_input_filename")
        self.file_service.set_file_id(input_file, "input")

    def onClose_x_cp2k_section_global_settings(self, backend, gIndex, section):
        # If the input file is available, parse it
        filepath = self.file_service.get_file_by_id("input")
        if filepath is not None:
            input_parser = CP2KInputParser(self.parser_context)
            input_parser.parse(filepath)
        else:
            logger.warning("The input file of the calculation could not be found.")

    def onClose_section_system(self, backend, gIndex, section):
        """Stores the index of the section method. Should always be 0, but
        let's get it dynamically just in case there's something wrong.
        """
        self.section_system_index = gIndex
        self.cache_service.addValue("number_of_atoms")
        self.cache_service.addArrayValues("configuration_periodic_dimensions")
        self.cache_service.addArrayValues("atom_labels")

    def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
        # Write the references to section_method and section_system
        backend.addValue('single_configuration_to_calculation_method_ref', self.section_method_index)
        backend.addValue('single_configuration_calculation_to_system_ref', self.section_system_index)

        scc_basis_id = backend.openSection("section_basis_set")

        # Basis kind
        self.cache_service.addValue("basis_set_kind")

        # Basis name
        basis_name = self.cache_service["basis_set_name"]
        if basis_name is not None:
            cutoff = self.cache_service["basis_set_planewave_cutoff"]
            if cutoff is not None:
                basis_name += "_PW_{}".format(cutoff)
                backend.addValue("basis_set_name", basis_name)

        # Gaussian mapping
        map_index_to_basis = []
        map_kind_number_to_basis_ref = self.cache_service["map_kind_number_to_basis_ref"]
        map_index_to_kind = self.cache_service["map_index_to_kind"]
        # print(map_kind_number_to_basis_ref)
        # print(map_index_to_kind)

        if map_index_to_kind is not None and map_kind_number_to_basis_ref is not None:
            indices = map_index_to_kind.keys()
            for index in sorted(indices):
                kind = map_index_to_kind[index]
                basis_ref = map_kind_number_to_basis_ref[kind]
                map_index_to_basis.append(basis_ref)
            map_index_to_basis = np.array(map_index_to_basis)
            backend.addArrayValues("mapping_section_basis_set_atom_centered", map_index_to_basis)

        # Cell dependent basis mapping
        self.cache_service.addValue("mapping_section_basis_set_cell_dependent")

        # If a cube file was found, get the filepath. Currently no information
        # is read from this file because there is no suitable metainfo for it.
        filename = self.cache_service["electron_density_filename"]
        if filename is not None:
            self.parser_context.file_service.get_absolute_path_to_file(filename)

        backend.closeSection("section_basis_set", scc_basis_id)

    #===========================================================================
    # adHoc functions
    def adHoc_x_cp2k_section_cell(self, parser):
        """Used to extract the cell information.
        """
        # Read the lines containing the cell vectors
        a_line = parser.fIn.readline()
        b_line = parser.fIn.readline()
        c_line = parser.fIn.readline()

        # Define the regex that extracts the components and apply it to the lines
        regex_string = r" CELL\| Vector \w \[angstrom\]:\s+({0})\s+({0})\s+({0})".format(self.regexs.float)
        regex_compiled = re.compile(regex_string)
        a_result = regex_compiled.match(a_line)
        b_result = regex_compiled.match(b_line)
        c_result = regex_compiled.match(c_line)

        # Convert the string results into a 3x3 numpy array
        cell = np.zeros((3, 3))
        cell[0, :] = [float(x) for x in a_result.groups()]
        cell[1, :] = [float(x) for x in b_result.groups()]
        cell[2, :] = [float(x) for x in c_result.groups()]

        # Push the results to cache
        self.cache_service["simulation_cell"] = cell
        self.cache_service["lattice_vectors"] = cell

    def adHoc_atom_forces(self, parser):
        """Used to extract the final atomic forces printed at the end of a
        calculation.
        """
        end_str = " SUM OF ATOMIC FORCES"
        end = False
        force_array = []

        # Loop through coordinates until the sum of forces is read
        while not end:
            line = parser.fIn.readline()
            if line.startswith(end_str):
                end = True
            else:
                forces = line.split()[-3:]
                forces = [float(x) for x in forces]
                force_array.append(forces)
        force_array = np.array(force_array)

        # If anything found, push the results to the correct section
        if len(force_array) != 0:
            # self.cache_service["atom_forces"] = force_array
            self.backend.addArrayValues("x_cp2k_atom_forces", force_array, unit="forceAu")

    def adHoc_stress_tensor(self, parser):
        """Used to extract the stress tensor printed at the end of a
        calculation.
        """
        row1 = [float(x) for x in parser.fIn.readline().split()[-3:]]
        row2 = [float(x) for x in parser.fIn.readline().split()[-3:]]
        row3 = [float(x) for x in parser.fIn.readline().split()[-3:]]
        stress_array = np.array([row1, row2, row3])
        parser.backend.addArrayValues("x_cp2k_stress_tensor", stress_array, unit="GPa")

    def adHoc_stress_calculation(self, parser):
        """Used to skip over the stress tensor calculation details.
        """
        end_line = " **************************** NUMERICAL STRESS END *****************************\n"
        finished = False
        while not finished:
            line = parser.fIn.readline()
            if line == end_line:
                finished = True

    def adHoc_stress_tensor_eigenpairs(self, parser):
        """Parses the stress tensor eigenpairs.
        """
        parser.fIn.readline()
        eigenvalues = np.array([float(x) for x in parser.fIn.readline().split()])
        parser.fIn.readline()
        row1 = [float(x) for x in parser.fIn.readline().split()]
        row2 = [float(x) for x in parser.fIn.readline().split()]
        row3 = [float(x) for x in parser.fIn.readline().split()]
        eigenvectors = np.array([row1, row2, row3])
        parser.backend.addArrayValues("x_cp2k_stress_tensor_eigenvalues", eigenvalues, unit="GPa")
        parser.backend.addArrayValues("x_cp2k_stress_tensor_eigenvectors", eigenvectors)

    def adHoc_single_point_converged(self, parser):
        """Called when the SCF cycle of a single point calculation has converged.
        """
        parser.backend.addValue("x_cp2k_quickstep_converged", True)

    def adHoc_single_point_not_converged(self, parser):
        """Called when the SCF cycle of a single point calculation did not converge.
        """
        parser.backend.addValue("x_cp2k_quickstep_converged", False)

    def adHoc_x_cp2k_section_quickstep_atom_information(self):
        """Used to extract the initial atomic coordinates and names in the
        Quickstep module.
        """
        def wrapper(parser):

            # Define the regex that extracts the information
            regex_string = r"\s+(\d+)\s+(\d+)\s+(\w+)\s+(\d+)\s+({0})\s+({0})\s+({0})".format(self.regexs.float)
            regex_compiled = re.compile(regex_string)

            match = True
            coordinates = []
            labels = []

            # Currently these three lines are not processed
            parser.fIn.readline()
            parser.fIn.readline()
            parser.fIn.readline()
            map_index_to_kind = {}

            while match:
                line = parser.fIn.readline()
                result = regex_compiled.match(line)

                if result:
                    match = True
                    kind_number = int(result.groups()[1])
                    element_number = int(result.groups()[3])
                    index = int(result.groups()[0])
                    map_index_to_kind[index] = kind_number
                    info = self.atom_kind_info[kind_number]
                    label = info["label"]
                    info["element_number"] = element_number
                    self.basis_set_info[kind_number]["element_number"] = element_number
                    labels.append(label)
                    coordinate = [float(x) for x in result.groups()[4:]]
                    coordinates.append(coordinate)
                else:
                    match = False
            coordinates = np.array(coordinates)
            labels = np.array(labels)
            self.cache_service["map_index_to_kind"] = map_index_to_kind

            # If anything found, push the results to the correct section
            if len(coordinates) != 0:
                self.cache_service["atom_positions"] = coordinates
                self.cache_service["atom_labels"] = labels

        return wrapper

    def adHoc_run_dir(self, metaname):
        def wrapper(parser):
            end_str = "\n"
            end = False
            path_array = []

            # Loop through lines until empty line is encountered
            while not end:
                line = parser.fIn.readline()
                if line == end_str or len(line) == 0:
                    end = True
                else:
                    path_part = line.split()[-1]
                    path_array.append(path_part)

            # Form the final path and push to backend
            path = "".join(path_array)
            parser.backend.addValue(metaname, path)

        return wrapper

    def adHoc_dft_plus_u(self, parser):
        self.test_electronic_structure_method = "DFT+U"

    def adHoc_mp2(self, parser):
        self.test_electronic_structure_method = "MP2"

    def adHoc_rpa(self, parser):
        self.test_electronic_structure_method = "RPA"

    def adHoc_print(self, msg):
        def wrapper(parser, groups):
            print(msg)
        return wrapper

    #=======================================================================
    # StarReActions
    def startReAction_save_cube_filename(self, backend, groups):
        filename = groups[0]
        self.cache_service["electron_density_filename"] = filename
