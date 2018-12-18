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

from __future__ import print_function
from __future__ import absolute_import
from builtins import next
from builtins import range
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.baseclasses import MainHierarchicalParser
import nomadcore.configurationreading
import nomadcore.csvparsing
from .commonparser import CP2KCommonParser
from nomadcore.caching_backend import CachingLevel
import logging
logger = logging.getLogger("nomad")


class CP2KGeoOptParser(MainHierarchicalParser):
    """Used to parse the CP2K calculation with run types:
        -GEO_OPT/GEOMETRY_OPTIMIZATION
    """
    def __init__(self, parser_context):
        """
        """
        super(CP2KGeoOptParser, self).__init__(parser_context)
        self.setup_common_matcher(CP2KCommonParser(parser_context))
        self.traj_iterator = None
        self.energy_reeval_quickstep = None

        #=======================================================================
        # Globally cached values
        self.cache_service.add("number_of_frames_in_sequence", 0)
        self.cache_service.add("frame_sequence_potential_energy", [])
        self.cache_service.add("frame_sequence_local_frames_ref", [])
        self.cache_service.add("geometry_optimization_method")

        #=======================================================================
        # Cache levels
        self.caching_levels.update({
            'x_cp2k_section_geometry_optimization_step': CachingLevel.ForwardAndCache,
            'x_cp2k_section_quickstep_calculation': CachingLevel.ForwardAndCache,
            'x_cp2k_section_geometry_optimization': CachingLevel.ForwardAndCache,
            # 'x_cp2k_section_geometry_optimization_energy_reevaluation': CachingLevel.ForwardAndCache,
        })

        #=======================================================================
        # SimpleMatchers
        self.geo_opt = SM(
            " ***                     STARTING GEOMETRY OPTIMIZATION                      ***".replace("*", "\*"),
            sections=["section_frame_sequence", "x_cp2k_section_geometry_optimization"],
            subMatchers=[
                SM( " ***                           CONJUGATE GRADIENTS                           ***".replace("*", "\*"),
                    adHoc=self.adHoc_conjugate_gradient(),
                    otherMetaInfo=["geometry_optimization_method"],
                ),
                SM( " ***                                   BFGS                                  ***".replace("*", "\*"),
                    adHoc=self.adHoc_bfgs(),
                    otherMetaInfo=["geometry_optimization_method"],
                ),
                SM( " ***                                 L-BFGS                                  ***".replace("*", "\*"),
                    adHoc=self.adHoc_bfgs(),
                    otherMetaInfo=["geometry_optimization_method"],
                ),
                # SM( "",
                    # forwardMatch=True,
                    # sections=["section_single_configuration_calculation", "section_system", "x_cp2k_section_geometry_optimization_step"],
                    # subMatchers=[
                        # self.cm.quickstep_calculation(),
                        # SM( " --------  Informations at step"),
                        # SM( "  Optimization Method        =\s+(?P<x_cp2k_optimization_method>{})".format(self.regexs.word)),
                        # SM( "  Total Energy               =\s+(?P<x_cp2k_optimization_energy__hartree>{})".format(self.regexs.float),
                            # otherMetaInfo=["frame_sequence_potential_energy"]
                        # ),
                    # ],
                    # otherMetaInfo=["atom_positions"],
                    # adHoc=self.adHoc_step(),
                # ),
                SM( " OPTIMIZATION STEP:",
                    endReStr="  Conv. in RMS gradients     =",
                    name="geooptstep",
                    repeats=True,
                    subMatchers=[
                        SM( "",
                            forwardMatch=True,
                            sections=["x_cp2k_section_geometry_optimization_step"],
                            otherMetaInfo=[
                                "atom_positions",
                            ],
                            subMatchers=[
                                # SM( "",
                                    # forwardMatch=True,
                                    # endReStr=" ***                 MNBRACK - NUMBER OF ENERGY EVALUATIONS :\s+{}\s+***".replace("*", "\*").format(self.regexs.int),
                                    # subMatchers=[
                                        # SM(" SCF WAVEFUNCTION OPTIMIZATION",
                                            # forwardMatch=True,
                                            # repeats=True,
                                            # subMatchers=[
                                                # self.cm.quickstep_calculation(),
                                            # ]
                                        # )
                                    # ]
                                # ),
                                # SM( "",
                                    # forwardMatch=True,
                                    # endReStr=" ***                 BRENT   - NUMBER OF ENERGY EVALUATIONS :\s+{}\s+***".replace("*", "\*").format(self.regexs.int),
                                    # subMatchers=[
                                        # SM(" SCF WAVEFUNCTION OPTIMIZATION",
                                            # forwardMatch=True,
                                            # repeats=True,
                                            # subMatchers=[
                                                # self.cm.quickstep_calculation(),
                                            # ]
                                        # )
                                    # ]
                                # ),
                                SM( " --------  Informations at step"),
                                SM( "  Optimization Method        =\s+(?P<x_cp2k_optimization_method>{})".format(self.regexs.word)),
                                SM( "  Total Energy               =\s+(?P<x_cp2k_optimization_energy__hartree>{})".format(self.regexs.float),
                                    otherMetaInfo=["frame_sequence_potential_energy"]
                                ),
                                SM( "  Real energy change         =\s+(?P<x_cp2k_optimization_energy_change__hartree>{})".format(self.regexs.float)),
                                SM( "  Decrease in energy         =\s+(?P<x_cp2k_optimization_energy_decrease>{})".format(self.regexs.word)),
                                SM( "  Used time                  =\s+(?P<x_cp2k_optimization_used_time>{})".format(self.regexs.float)),
                                SM( "  Max. step size             =\s+(?P<x_cp2k_optimization_max_step_size__bohr>{})".format(self.regexs.float)),
                                SM( "  Conv. limit for step size  =\s+(?P<x_cp2k_optimization_step_size_convergence_limit__bohr>{})".format(self.regexs.float),
                                    otherMetaInfo=["geometry_optimization_geometry_change"]
                                ),
                                SM( "  Convergence in step size   =\s+(?P<x_cp2k_optimization_step_size_convergence>{})".format(self.regexs.word)),
                                SM( "  RMS step size              =\s+(?P<x_cp2k_optimization_rms_step_size__bohr>{})".format(self.regexs.float)),
                                SM( "  Convergence in RMS step    =\s+(?P<x_cp2k_optimization_rms_step_size_convergence>{})".format(self.regexs.word)),
                                SM( "  Max. gradient              =\s+(?P<x_cp2k_optimization_max_gradient__bohr_1hartree>{})".format(self.regexs.float)),
                                SM( "  Conv. limit for gradients  =\s+(?P<x_cp2k_optimization_gradient_convergence_limit__bohr_1hartree>{})".format(self.regexs.float),
                                    otherMetaInfo=["geometry_optimization_threshold_force"]
                                ),
                                SM( "  Conv. for gradients        =\s+(?P<x_cp2k_optimization_max_gradient_convergence>{})".format(self.regexs.word)),
                                SM( "  RMS gradient               =\s+(?P<x_cp2k_optimization_rms_gradient__bohr_1hartree>{})".format(self.regexs.float)),
                                SM( "  Conv. in RMS gradients     =\s+(?P<x_cp2k_optimization_rms_gradient_convergence>{})".format(self.regexs.word)),
                            ],
                            # adHoc=self.adHoc_step()
                        ),
                    ]
                ),
                SM( " ***                    GEOMETRY OPTIMIZATION COMPLETED                      ***".replace("*", "\*"),
                    adHoc=self.adHoc_geo_opt_converged(),
                    otherMetaInfo=["geometry_optimization_converged"]
                ),
                SM( "                    Reevaluating energy at the minimum",
                    # sections=["x_cp2k_section_geometry_optimization_energy_reevaluation"],
                    subMatchers=[
                        self.cm.quickstep_calculation(),
                        # SM("",
                            # adHoc=self.adHoc_save_energy_reeval_quickstep()
                        # )
                    ],
                    # adHoc=self.adHoc_save_energy_reeval_quickstep()
                ),
                # SM( "",
                    # forwardMatch=True,
                    # adHoc=self.adHoc_save_energy_reeval_quickstep()
                # )
            ],
        )

        # Compose root matcher according to the run type. This way the
        # unnecessary regex parsers will not be compiled and searched. Saves
        # computational time.
        self.root_matcher = SM("",
            forwardMatch=True,
            sections=["section_run", "section_sampling_method"],
            subMatchers=[
                SM( "",
                    forwardMatch=True,
                    sections=["section_method"],
                    subMatchers=[
                        self.cm.header(),
                        self.cm.quickstep_header(),
                    ],
                ),
                self.geo_opt,
                self.cm.footer(),
            ]
        )

    #===========================================================================
    # onClose triggers
    def onClose_x_cp2k_section_geometry_optimization(self, backend, gIndex, section):

        # Get the re-evaluated energy and add it to frame_sequence_potential_energy
        reeval_quickstep = self.energy_reeval_quickstep
        if reeval_quickstep is not None:
            energy = reeval_quickstep.get_latest_value("x_cp2k_energy_total")
            if energy is not None:
                self.cache_service["frame_sequence_potential_energy"].append(energy)

        # Push values from cache
        self.cache_service.addArrayValues("frame_sequence_potential_energy")
        self.cache_service.addValue("geometry_optimization_method")
        self.backend.addValue("frame_sequence_to_sampling_ref", 0)

        # Get the optimization convergence criteria from the last optimization
        # step
        section.add_latest_value([
            "x_cp2k_section_geometry_optimization_step",
            "x_cp2k_optimization_step_size_convergence_limit"],
            "geometry_optimization_geometry_change",
        )
        section.add_latest_value([
            "x_cp2k_section_geometry_optimization_step",
            "x_cp2k_optimization_gradient_convergence_limit"],
            "geometry_optimization_threshold_force",
        )

        # Push the information into single configuration and system
        steps = section["x_cp2k_section_geometry_optimization_step"]
        each = self.cache_service["each_geo_opt"]
        add_last = False
        add_last_setting = self.cache_service["traj_add_last"]
        if add_last_setting == "NUMERIC" or add_last_setting == "SYMBOLIC":
            add_last = True

        # Push the trajectory
        n_steps = len(steps) + 1
        last_step = n_steps - 1

        # First push the original system geometry
        # print(self.cache_service["atom_positions"])
        for i_step in range(n_steps):
            singleId = backend.openSection("section_single_configuration_calculation")
            systemId = backend.openSection("section_system")

            if self.traj_iterator is not None:
                if (i_step + 1) % each == 0 or (i_step == last_step and add_last):
                    try:
                        pos = next(self.traj_iterator)
                    except StopIteration:
                        logger.error(
                            "Could not get next geometry from an external"
                            " file. It seems that the number of optimization "
                            "steps in the CP2K output file doesn't match the "
                            "number of steps found in the external trajectory "
                            "file."
                        )
                    else:
                        backend.addArrayValues("atom_positions", pos, unit="angstrom")
            backend.closeSection("section_system", systemId)
            backend.closeSection("section_single_configuration_calculation", singleId)

        self.cache_service.addArrayValues("frame_sequence_local_frames_ref")
        backend.addValue("number_of_frames_in_sequence", n_steps)

    def onClose_section_sampling_method(self, backend, gIndex, section):
        self.backend.addValue("sampling_method", "geometry_optimization")

    def onClose_x_cp2k_section_quickstep_calculation(self, backend, gIndex, section):
        self.energy_reeval_quickstep = section

    def onClose_x_cp2k_section_geometry_optimization_step(self, backend, gIndex, section):
        energy = section["x_cp2k_optimization_energy"]
        if energy is not None:
            self.cache_service["frame_sequence_potential_energy"].append(energy[0])

    def onClose_section_system(self, backend, gIndex, section):
        self.cache_service.addArrayValues("simulation_cell", unit="angstrom")
        self.cache_service.addArrayValues("lattice_vectors", unit="angstrom")

    def onClose_section_method(self, backend, gIndex, section):
        traj_file = self.file_service.get_file_by_id("trajectory")
        traj_format = self.cache_service["trajectory_format"]
        if traj_format is not None and traj_file is not None:

            # Use special parsing for CP2K pdb files because they don't follow the proper syntax
            if traj_format == "PDB":
                self.traj_iterator = nomadcore.csvparsing.iread(traj_file, columns=[3, 4, 5], start="CRYST", end="END")
            else:
                try:
                    self.traj_iterator = nomadcore.configurationreading.iread(traj_file)
                except ValueError:
                    pass

    def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
        self.cache_service["frame_sequence_local_frames_ref"].append(gIndex)

    #===========================================================================
    # adHoc functions
    def adHoc_geo_opt_converged(self):
        """Called when the geometry optimization converged.
        """
        def wrapper(parser):
            parser.backend.addValue("geometry_optimization_converged", True)
        return wrapper

    def adHoc_geo_opt_not_converged(self):
        """Called when the geometry optimization did not converge.
        """
        def wrapper(parser):
            parser.backend.addValue("geometry_optimization_converged", False)
        return wrapper

    def adHoc_conjugate_gradient(self):
        """Called when conjugate gradient method is used.
        """
        def wrapper(parser):
            self.cache_service["geometry_optimization_method"] = "conjugate_gradient"
        return wrapper

    def adHoc_bfgs(self):
        """Called when conjugate gradient method is used.
        """
        def wrapper(parser):
            self.cache_service["geometry_optimization_method"] = "bfgs"
        return wrapper

    # def adHoc_save_energy_reeval_quickstep(self):
        # def wrapper(parser):
            # section_managers = parser.backend.sectionManagers
            # section_run_manager = section_managers["section_run"]
            # section_run = section_run_manager.openSections[0]
            # print section_run.subSectionValues
            # # quickstep = section_run.get_latest_value("x_cp2k_section_quickstep_calculation")
            # # print quickstep
            # # self.energy_reeval_quickstep = quickstep
        # return wrapper

    def debug(self):
        def wrapper(parser):
            print("DEBUG")
        return wrapper
