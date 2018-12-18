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
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.baseclasses import MainHierarchicalParser
from cp2kparser.versions.cp2k262.singlepointforceparser import CP2KSinglePointForceParser
from nomadcore.caching_backend import CachingLevel
from cp2kparser.versions.cp2k262.commonparser import CP2KCommonParser
import logging
logger = logging.getLogger("nomad")


class CP2KSinglePointParser(MainHierarchicalParser):
    """The main parser class. Used to parse the CP2K calculation with run types:
        -ENERGY
        -ENERGY_FORCE
    """
    def __init__(self, parser_context):
        """
        """
        super(CP2KSinglePointParser, self).__init__(parser_context)
        self.setup_common_matcher(CP2KCommonParser(parser_context))

        #=======================================================================
        # Cache levels
        self.caching_levels.update({
            'x_cp2k_energy_total_scf_iteration': CachingLevel.ForwardAndCache,
            'x_cp2k_energy_XC_scf_iteration': CachingLevel.ForwardAndCache,
            'x_cp2k_energy_change_scf_iteration': CachingLevel.ForwardAndCache,
            'x_cp2k_stress_tensor': CachingLevel.ForwardAndCache,
            'x_cp2k_section_stress_tensor': CachingLevel.ForwardAndCache,
        })

        #=======================================================================
        # SimpleMatchers
        self.root_matcher = SM("",
            forwardMatch=True,
            sections=['section_run', "section_single_configuration_calculation", "section_system", "section_method"],
            otherMetaInfo=["atom_forces"],
            subMatchers=[
                self.cm.header(),
                self.cm.quickstep_header(),
                self.cm.quickstep_calculation(),
                self.cm.footer(),
            ]
        )

    #===========================================================================
    # onClose triggers
    def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
        """
        """
        # If the force file for a single point calculation is available, and
        # the forces were not parsed fro the output file, parse the separate
        # file
        if section["atom_forces"] is None:
            force_file = self.file_service.get_file_by_id("force_file_single_point")
            if force_file is not None:
                force_parser = CP2KSinglePointForceParser(self.parser_context)
                force_parser.parse(force_file)
            else:
                logger.warning("The file containing the forces printed by ENERGY_FORCE calculation could not be found.")

        # Only in the single configuration calculations the number of scf
        # iterations is given. E.g. in geometry optimization there are multiple
        # scf calculations so this loses it's meaning sort of.
        self.cache_service.addValue("number_of_scf_iterations")

    def onClose_x_cp2k_section_scf_iteration(self, backend, gIndex, section):
        """Keep track of how many SCF iteration are made."""
        self.cache_service["number_of_scf_iterations"] += 1
        gId = backend.openSection("section_scf_iteration")
        section.add_latest_value("x_cp2k_energy_total_scf_iteration", "energy_total_scf_iteration")
        section.add_latest_value("x_cp2k_energy_XC_scf_iteration", "energy_XC_scf_iteration")
        section.add_latest_value("x_cp2k_energy_change_scf_iteration", "energy_change_scf_iteration")
        backend.closeSection("section_scf_iteration", gId)

    def onClose_x_cp2k_section_quickstep_calculation(self, backend, gIndex, section):
        """"""
        section.add_latest_value("x_cp2k_energy_total", "energy_total")
        section.add_latest_value("x_cp2k_electronic_kinetic_energy", "electronic_kinetic_energy")
        section.add_latest_value("x_cp2k_quickstep_converged", "single_configuration_calculation_converged")
        section.add_latest_array_values("x_cp2k_atom_forces", "atom_forces")

    def onClose_x_cp2k_section_stress_tensor(self, backend, gIndex, section):
        """"""
        gId = backend.openSection("section_stress_tensor")
        section.add_latest_array_values("x_cp2k_stress_tensor", "stress_tensor")
        backend.closeSection("section_stress_tensor", gId)

    def onClose_section_system(self, backend, gIndex, section):
        """Stores the index of the section method. Should always be 0, but
        let's get it dynamically just in case there's something wrong.
        """
        self.cache_service.addArrayValues("atom_positions", unit="angstrom")
        self.cache_service.addArrayValues("simulation_cell", unit="angstrom")
        self.cache_service.addArrayValues("lattice_vectors", unit="angstrom")

    #===========================================================================
    # adHoc functions
