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

from future import standard_library
standard_library.install_aliases()
from builtins import object
import os
import re
from io import StringIO
import logging
import pickle
import numpy as np
from nomadcore.baseclasses import AbstractBaseParser
import nomadcore.configurationreading
from cp2kparser.generic.inputparsing import metainfo_data_prefix, metainfo_section_prefix
from pint import UnitRegistry
import ase

ureg = UnitRegistry()
logger = logging.getLogger("nomad")


class CP2KInputParser(AbstractBaseParser):
    """Used to parse out a CP2K input file.

    CP2K offers a complete structure for the input in an XML file, which can be
    printed with the command cp2k --xml. This XML file has been preparsed into
    a native python object ('CP2KInput' class found in generic.inputparsing)
    and stored in a python pickle file. It e.g. contains all the default values
    that are often needed as they are used if the user hasn't specified a
    settings in the input. This XML file is used to get the default values
    because it is rather cumbersome to hard code them in the parser itself,
    especially if there will be lot's of them. Hard coded values will also be
    more error prone, and would have to be checked for each parser version.

    CP2K input supports including other input files and also
    supports variables. This is currently not supported, but may be added at
    some point.
    """
    def __init__(self, parser_context):
        """
        Attributes:
            input_tree: The input structure for this version of CP2K. The
                structure is already present, in this module it will be filled with
                data found from the input file.
            input_lines: List of preprocessed lines in the input. Here all the
                variables have been stated explicitly and the additional input files have
                been merged.
        """
        super(CP2KInputParser, self).__init__(parser_context)
        self.input_tree = None
        self.input_lines = None
        self.unit_mapping = {
            # Distance
            "BOHR": "bohr",
            "M": "m",
            "PM": "pm",
            "NM": "nm",
            "ANGSTROM": "angstrom",
            # Time
            "S": "s",
            "FS": "fs",
            "PS": "ps",
            "AU_T": "(hbar/hartree)",
            "WAVENUMBER_T": None,
        }

        #=======================================================================
        # Cached values
        self.cache_service.add("configuration_periodic_dimensions", single=False, update=False)
        self.cache_service.add("trajectory_format")
        self.cache_service.add("trajectory_unit")
        self.cache_service.add("velocity_format")
        self.cache_service.add("velocity_unit")
        self.cache_service.add("vel_add_last")
        self.cache_service.add("each_geo_opt")
        self.cache_service.add("traj_add_last")

    def parse(self, filepath):

        # Preprocess to spell out variables and to include stuff from other
        # files
        self.preprocess_input(filepath)

        # Gather the information from the input file
        self.fill_input_tree(filepath)

        # Parse everything in the input to cp2k specific metadata
        self.fill_metadata()

        # If the print level is low, parse the system geometry from the input
        # file. If keywords that manipulate the system are enabled, the
        # geometry is not parsed because the manipulations can be hard to
        # reproduce.
        print_level = self.input_tree.get_keyword("GLOBAL/PRINT_LEVEL")
        multiple_unit_cell = self.input_tree.get_keyword("FORCE_EVAL/SUBSYS/MULTIPLE_UNIT_CELL")
        top_multiple_unit_cell = self.input_tree.get_keyword("FORCE_EVAL/SUBSYS/TOPOLOGY/MULTIPLE_UNIT_CELL")
        centering = self.input_tree.get_section("FORCE_EVAL/SUBSYS/TOPOLOGY/CENTER_COORDINATES").section_parameter.get_value()
        center_point = self.input_tree.get_keyword("FORCE_EVAL/SUBSYS/TOPOLOGY/CENTER_COORDINATES/CENTER_POINT")
        if multiple_unit_cell is None:
            multiple_unit_cell = np.array([1, 1, 1])

        if print_level == "LOW" \
           and np.array_equal(multiple_unit_cell, [1, 1, 1]) \
           and np.array_equal(top_multiple_unit_cell, [1, 1, 1]):

            cell_section = self.input_tree.get_section("FORCE_EVAL/SUBSYS/CELL")
            a_obj = cell_section.get_keyword_object("A")
            b_obj = cell_section.get_keyword_object("B")
            c_obj = cell_section.get_keyword_object("C")
            a_val = a_obj.get_value()
            b_val = b_obj.get_value()
            c_val = c_obj.get_value()

            abc_obj = cell_section.get_keyword_object("ABC")
            angles_obj = cell_section.get_keyword_object("ALPHA_BETA_GAMMA")
            abc_val = abc_obj.get_value()
            angles_val = angles_obj.get_value()

            # Cell given as three vectors
            if a_val is not None and b_val is not None and c_val is not None:
                a_unit = a_obj.get_unit()
                b_unit = b_obj.get_unit()
                c_unit = c_obj.get_unit()
                a = (a_val * ureg(a_unit)).to("angstrom").magnitude
                b = (b_val * ureg(b_unit)).to("angstrom").magnitude
                c = (c_val * ureg(c_unit)).to("angstrom").magnitude
                cell = np.stack((a, b, c), axis=0)
                self.cache_service["simulation_cell"] = cell
                self.cache_service["lattice_vectors"] = cell
            # Cell given as lengths and angles
            elif abc_val is not None and angles_val is not None:
                abc_unit = abc_obj.get_unit()
                angles_unit = angles_obj.get_unit()
                abc = (abc_val * ureg(abc_unit)).to("angstrom").magnitude
                angles = (angles_val * ureg(angles_unit)).to("degree").magnitude
                parameters = np.hstack((abc, angles))
                cell = ase.geometry.cellpar_to_cell(parameters)
                self.cache_service["simulation_cell"] = cell
                self.cache_service["lattice_vectors"] = cell

            # Parse the coordinates from input file or xyz file
            coord = self.input_tree.get_section("FORCE_EVAL/SUBSYS/COORD")
            coord_filename = self.input_tree.get_keyword("FORCE_EVAL/SUBSYS/TOPOLOGY/COORD_FILE_NAME")
            coord_format = self.input_tree.get_keyword("FORCE_EVAL/SUBSYS/TOPOLOGY/COORD_FILE_FORMAT")
            scaled = coord.get_keyword("SCALED")
            coords = coord.default_keyword.value
            coord_unit = coord.get_keyword("UNIT")

            final_pos = None
            final_lab = None
            if len(coords) != 0:
                lines = coords.split("\n")
                pos = []
                labels = []
                for line in lines:
                    if len(line) != 0:
                        parts = line.split()
                        i_lab = parts[0]
                        i_x = float(parts[1])
                        i_y = float(parts[2])
                        i_z = float(parts[3])
                        labels.append(i_lab)
                        pos.append([i_x, i_y, i_z])

                final_pos = np.array(pos)
                final_lab = np.array(labels)

            elif coord_filename is not None:

                extension_map = {
                    "XYZ": "xyz",
                }
                coord_extension = extension_map.get(coord_format)
                if coord_extension is not None:

                    # We first remove the comment line from the XYZ file,
                    # because ASE cannot properly handle comments that have ASE
                    # specific data that is malformed.
                    abs_filename = self.file_service.get_absolute_path_to_file(coord_filename)
                    stringio = StringIO()
                    with open(abs_filename, "r") as fin:
                        for i_line, line in enumerate(fin):
                            if i_line == 1:
                                stringio.write("\n")
                            else:
                                stringio.write(line)
                    try:
                        atoms = ase.io.read(stringio, format="xyz")
                        final_pos = atoms.get_positions()
                        final_lab = atoms.get_chemical_symbols()
                    except StopIteration:
                        logger.error("Error in reading the structure file specified in the input file.")
                else:
                    logger.error(
                        "Could not read the structure file specified in the "
                        "input file because the file format is not supported."
                    )

            # Make the positions cartesian if they are scaled
            if final_pos is not None and final_lab is not None:
                if scaled == "T" or scaled == "TRUE":
                    atoms = ase.Atoms(
                        scaled_positions=final_pos,
                        symbols=final_lab,
                        cell=cell,
                    )
                    final_pos = atoms.get_positions()

                # Center the atoms if requested.
                if centering == "T" or centering == "TRUE":
                    atoms = ase.Atoms(
                        positions=final_pos,
                        symbols=final_lab,
                        cell=cell,
                    )
                    atoms.center(about=center_point)
                    final_pos = atoms.get_positions()

                coord_unit = self.get_pint_unit_string(coord_unit)
                final_pos = (final_pos * ureg(coord_unit)).to("angstrom").magnitude
                self.cache_service["atom_positions"] = final_pos
                self.cache_service["atom_labels"] = final_lab

        # Parse the used XC_functionals and their parameters
        xc = self.input_tree.get_section("FORCE_EVAL/DFT/XC/XC_FUNCTIONAL")
        if xc is not None:
            xc_list = []

            class XCFunctional(object):
                def __init__(self, name, weight=1, parameters=None):
                    self.name = name
                    self.weight = weight
                    self.parameters = parameters

            # First see if a functional has been specified in the section parameter
            section_parameter = xc.section_parameter.value
            if section_parameter is not None and section_parameter != "NO_SHORTCUT":

                if section_parameter == "BLYP":
                    xc_list.append(XCFunctional("GGA_X_B88"))
                    xc_list.append(XCFunctional("GGA_C_LYP"))

                elif section_parameter == "LDA" or section_parameter == "PADE":
                    xc_list.append(XCFunctional("LDA_XC_TETER93"))

                elif section_parameter == "PBE":
                    xc_list.append(XCFunctional("GGA_X_PBE"))
                    xc_list.append(XCFunctional("GGA_C_PBE"))

                elif section_parameter == "OLYP":
                    xc_list.append(XCFunctional("GGA_X_OPTX"))
                    xc_list.append(XCFunctional("GGA_C_LYP"))

                elif section_parameter == "HCTH120":
                    xc_list.append(XCFunctional("GGA_XC_HCTH_120"))

                elif section_parameter == "PBE0":
                    xc_list.append(XCFunctional("HYB_GGA_XC_PBEH"))

                elif section_parameter == "B3LYP":
                    xc_list.append(XCFunctional("HYB_GGA_XC_B3LYP"))

                elif section_parameter == "TPSS":
                    xc_list.append(XCFunctional("MGGA_X_TPSS"))
                    xc_list.append(XCFunctional("MGGA_C_TPSS"))

                else:
                    logger.warning("Unknown XC functional given in XC_FUNCTIONAL section parameter.")

            # Otherwise one has to look at the individual functional settings
            else:
                pbe = xc.get_subsection("PBE")
                if pbe is not None:
                    if pbe.accessed:
                        sp = pbe.get_section_parameter()
                        if sp == "T":
                            parametrization = pbe.get_keyword("PARAMETRIZATION", allow_default=True)
                            scale_x = pbe.get_keyword("SCALE_X", allow_default=True)
                            scale_c = pbe.get_keyword("SCALE_C", allow_default=True)
                            if parametrization == "ORIG":
                                xc_list.append(XCFunctional("GGA_X_PBE", scale_x))
                                xc_list.append(XCFunctional("GGA_C_PBE", scale_c))
                            elif parametrization == "PBESOL":
                                xc_list.append(XCFunctional("GGA_X_PBE_SOL", scale_x))
                                xc_list.append(XCFunctional("GGA_C_PBE_SOL", scale_c))
                            elif parametrization == "REVPBE":
                                xc_list.append(XCFunctional("GGA_X_PBE_R", scale_x))
                                xc_list.append(XCFunctional("GGA_C_PBE", scale_c))
                tpss = xc.get_subsection("TPSS")
                if tpss is not None:
                    if tpss.accessed:
                        sp = tpss.get_section_parameter()
                        if sp == "T":
                            scale_x = tpss.get_keyword("SCALE_X", allow_default=True)
                            scale_c = tpss.get_keyword("SCALE_C", allow_default=True)
                            xc_list.append(XCFunctional("MGGA_X_TPSS", scale_x))
                            xc_list.append(XCFunctional("MGGA_C_TPSS", scale_c))

            # Sort the functionals alphabetically by name
            xc_list.sort(key=lambda x: x.name)
            xc_summary = ""

            # For every defined functional, stream the information to the
            # backend and construct the summary string
            for i, functional in enumerate(xc_list):

                gId = self.backend.openSection("section_XC_functionals")
                self.backend.addValue("XC_functional_name", functional.name)
                self.backend.addValue("XC_functional_weight", functional.weight)
                if functional.parameters is not None:
                    pass
                self.backend.closeSection("section_XC_functionals", gId)

                if i != 0:
                    xc_summary += "+"
                xc_summary += "{}*{}".format(functional.weight, functional.name)
                if functional.parameters is not None:
                    xc_summary += ":{}".format()

            # Stream summary
            if xc_summary is not "":
                self.backend.addValue("XC_functional", xc_summary)

        #=======================================================================
        # Cell periodicity
        periodicity = self.input_tree.get_keyword("FORCE_EVAL/SUBSYS/CELL/PERIODIC")
        if periodicity is not None:
            periodicity = periodicity.upper()
            periodicity_list = ("X" in periodicity, "Y" in periodicity, "Z" in periodicity)
            self.cache_service["configuration_periodic_dimensions"] = np.asarray(periodicity_list)
        else:
            logger.warning("Could not determine cell periodicity from FORCE_EVAL/SUBSYS/CELL/PERIODIC")

        #=======================================================================
        # Single point force file name
        self.setup_force_file_name()

        #=======================================================================
        # Trajectory file name
        self.setup_trajectory_file_name()

        #=======================================================================
        # Trajectory file format
        self.cache_service["trajectory_format"] = self.input_tree.get_keyword("MOTION/PRINT/TRAJECTORY/FORMAT")
        self.cache_service["traj_add_last"] = self.input_tree.get_keyword("MOTION/PRINT/TRAJECTORY/ADD_LAST")
        traj_unit = self.input_tree.get_keyword("MOTION/PRINT/TRAJECTORY/UNIT")
        pint_traj_unit = self.get_pint_unit_string(traj_unit)
        self.cache_service["trajectory_unit"] = pint_traj_unit

        #=======================================================================
        # Velocity file format
        self.cache_service["velocity_format"] = self.input_tree.get_keyword("MOTION/PRINT/VELOCITIES/FORMAT")
        self.cache_service["vel_add_last"] = self.input_tree.get_keyword("MOTION/PRINT/VELOCITIES/ADD_LAST")
        vel_unit = self.input_tree.get_keyword("MOTION/PRINT/VELOCITIES/UNIT")
        pint_vel_unit = self.get_pint_unit_string(vel_unit)
        self.cache_service["velocity_unit"] = pint_vel_unit

        #=======================================================================
        # OLD: This information is retrieved from output file
        # See if some more exotic calculation is requested (e.g. MP2, DFT+U, GW, RPA)
        # Search for a WF_CORRELATION section
        # correlation = self.input_tree.get_section("FORCE_EVAL/DFT/XC/WF_CORRELATION")
        # method = "DFT"
        # if correlation.accessed:
            # method = correlation.get_keyword_value_raw("METHOD")
            # if method != "NONE":
                # # Can't really decide which method used (MP2, RPA, GW)
                # method = None

        # # Search for DFT+U settings
        # kinds = self.input_tree.get_section_list("FORCE_EVAL/SUBSYS/KIND")
        # for kind in kinds:
            # dft_u = kind.get_subsection("DFT_PLUS_U")
            # if dft_u.accessed:
                # method = "DFT+U"

        # self.cache_service["electronic_structure_method"] = method

        #=======================================================================
        # Stress tensor calculation method
        stress_tensor_method = self.input_tree.get_keyword("FORCE_EVAL/STRESS_TENSOR")
        if stress_tensor_method != "NONE":
            mapping = {
                "NUMERICAL": "Numerical",
                "ANALYTICAL": "Analytical",
                "DIAGONAL_ANALYTICAL": "Diagonal analytical",
                "DIAGONAL_NUMERICAL": "Diagonal numerical",
            }
            stress_tensor_method = mapping.get(stress_tensor_method)
            if stress_tensor_method is not None:
                self.backend.addValue("stress_tensor_method", stress_tensor_method)

    def normalize_x_cp2k_path(self, path):
        """The paths in CP2K input can be given in many ways. This function
        tries to normalize these forms into a valid path.
        """
        # Path is exactly as given
        if path.startswith("="):
            normalized_path = path[1:]
        # Path is relative, no project name added
        elif re.match(r"./", path):
            normalized_path = path
        # Path is relative, project name added
        else:
            project_name = self.input_tree.get_keyword("GLOBAL/PROJECT_NAME")
            if path:
                normalized_path = "{}-{}".format(project_name, path)
            else:
                normalized_path = project_name
        return normalized_path

    def get_pint_unit_string(self, cp2k_unit_string):
        """Translate the CP2K unit definition into a valid pint unit.
        """
        units = re.split('[\^\-\+\*\d]+', cp2k_unit_string)
        for unit in units:
            if unit == "":
                continue
            pint_unit = self.unit_mapping.get(unit.upper())
            if pint_unit is None:
                return None
            cp2k_unit_string = cp2k_unit_string.replace(unit, pint_unit)
        return cp2k_unit_string

    def setup_force_file_name(self):
        """Setup the force file path.
        """
        force_file = self.input_tree.get_keyword("FORCE_EVAL/PRINT/FORCES/FILENAME")
        extension = "xyz"
        if force_file is not None and force_file != "__STD_OUT__":
            normalized_path = self.normalize_x_cp2k_path(force_file)
            final_path = "{}-1_0.{}".format(normalized_path, extension)
            self.file_service.set_file_id(final_path, "force_file_single_point")

    def setup_trajectory_file_name(self):
        """Setup the trajectory file path.
        """
        traj_format = self.input_tree.get_keyword("MOTION/PRINT/TRAJECTORY/FORMAT")
        traj_filename = self.input_tree.get_keyword("MOTION/PRINT/TRAJECTORY/FILENAME")
        self.cache_service["each_geo_opt"] = self.input_tree.get_keyword("MOTION/PRINT/TRAJECTORY/EACH/GEO_OPT")
        if traj_filename is None:
            traj_filename = ""
        extension_map = {
            "XYZ": "xyz",
            "XMOL": "xyz",
            "ATOMIC": "xyz",
            "PDB": "pdb",
            "DCD": "dcd",
        }
        extension = extension_map.get(traj_format)
        if extension is None:
            logger.error("Unknown file format '{}' for CP2K trajectory file ".format(traj_format))
            return
        normalized_path = self.normalize_x_cp2k_path(traj_filename)
        final_path = "{}-pos-1.{}".format(normalized_path, extension)
        self.file_service.set_file_id(final_path, "trajectory")

    def fill_input_tree(self, file_path):
        """Parses a CP2K input file into an object tree.

        Return an object tree represenation of the input augmented with the
        default values and lone keyword values from the x_cp2k_input.xml file
        which is version specific. Keyword aliases are also mapped to the same
        data.

        The cp2k input is largely case-insensitive. In the input tree, we wan't
        only one standard way to name things, so all section names and section
        parameters will be transformed into upper case.

        To query the returned tree use the following functions:
            get_keyword("GLOBAL/PROJECT_NAME")
            get_parameter("GLOBAL/PRINT")
            get_default_keyword("FORCE_EVAL/SUBSYS/COORD")

        Args:
            : A string containing the contents of a CP2K input file. The
            input file can be stored as string as it isn't that big.

        Returns:
            The input as an object tree.
        """

        self.setup_version(self.parser_context.version_id)
        section_stack = []
        self.input_tree.root_section.accessed = True

        for line in self.input_lines:

            # Remove comments and whitespaces
            line = line.split('!', 1)[0].split('#', 1)[0].strip()

            # Skip empty lines
            if len(line) == 0:
                continue

            # Section ends
            if line.upper().startswith('&END'):
                section_stack.pop()
            # Section starts
            elif line[0] == '&':
                parts = line.split(' ', 1)
                name = parts[0][1:].upper()
                section_stack.append(name)

                # Form the path
                path = ""
                for index, item in enumerate(section_stack):
                    if index != 0:
                        path += '/'
                    path += item

                # Mark the section as accessed.
                self.input_tree.set_section_accessed(path)

                # Save the section parameters
                if len(parts) > 1:
                    self.input_tree.set_parameter(path, parts[1].strip().upper())

            # Ignore variables and includes that might still be here for some
            # reason
            elif line.upper().startswith('@'):
                continue

            # Contents (keywords, default keywords)
            else:
                split = line.split(None, 1)
                if len(split) <= 1:
                    keyword_value = ""
                else:
                    keyword_value = split[1]
                keyword_name = split[0].upper()
                self.input_tree.set_keyword(path + "/" + keyword_name, keyword_value, line)

    def fill_metadata(self):
        """Goes through the input data and pushes everything to the
        backend.
        """
        name_stack = []
        self.fill_metadata_recursively(self.input_tree.root_section, name_stack)

    def fill_metadata_recursively(self, section, name_stack):
        """Recursively goes through the input sections and pushes everything to the
        backend.
        """
        if not section.accessed:
            return

        if section.name == "CP2K_INPUT":
            path = "x_cp2k_section_input"
        else:
            name_stack.append(section.name)
            path = metainfo_section_prefix + "{}".format(".".join(name_stack))

        not_section_path = metainfo_data_prefix + "{}".format(".".join(name_stack))

        gid = self.backend.openSection(path)

        # Keywords
        for default_name in section.default_keyword_names:
            keywords = section.keywords.get(default_name)
            for keyword in keywords:
                if keyword.value is not None:
                    name = "{}.{}".format(not_section_path, keyword.default_name)
                    self.backend.addValue(name, keyword.value)

        # Section parameter
        section_parameter = section.section_parameter
        if section_parameter is not None:
            name = "{}.SECTION_PARAMETERS".format(not_section_path)
            if section_parameter.value is not None:
                self.backend.addValue(name, section_parameter.value)

        # Default keyword
        default_keyword = section.default_keyword
        if default_keyword is not None:

            name = "{}.DEFAULT_KEYWORD".format(not_section_path)
            self.backend.addValue(name, default_keyword.value)

        # Subsections
        for name, subsections in section.sections.items():
            for subsection in subsections:
                self.fill_metadata_recursively(subsection, name_stack)

        self.backend.closeSection(path, gid)

        if section.name != "CP2K_INPUT":
            name_stack.pop()

    def setup_version(self, version_number):
        """ The pickle file which contains preparsed data from the
        x_cp2k_input.xml is version specific. By calling this function before
        parsing the correct file can be found.
        """
        pickle_path = os.path.dirname(__file__) + "/input_data/cp2k_input_tree.pickle"
        input_tree_pickle_file = open(pickle_path, 'rb')
        self.input_tree = pickle.load(input_tree_pickle_file)

    def preprocess_input(self, filepath):
        """Preprocess the input file. Concatenate .inc files into the main
        input file and explicitly state all variables.
        """
        # Read the input file into memory. It shouldn't be that big so we can
        # do this easily
        input_lines = []
        with open(filepath, "r") as f:
            for line in f:
                input_lines.append(line.strip())

        # Merge include files to input
        extended_input = input_lines[:]  # Make a copy
        i_line = 0
        for line in input_lines:
            if line.startswith("@INCLUDE") or line.startswith("@include"):
                split = line.split(None, 1)
                includepath = split[1]
                basedir = os.path.dirname(filepath)
                filepath = os.path.join(basedir, includepath)
                filepath = os.path.abspath(filepath)
                if not os.path.isfile(filepath):
                    logger.warning("Could not find the include file '{}' stated in the CP2K input file. Continuing without it.".format(filepath))
                    continue

                # Get the content from include file
                included_lines = []
                with open(filepath, "r") as includef:
                    for line in includef:
                        included_lines.append(line.strip())
                    del extended_input[i_line]
                    extended_input[i_line:i_line] = included_lines
                    i_line += len(included_lines)
            i_line += 1

        # Gather the variable definitions
        variables = {}
        input_set_removed = []
        for i_line, line in enumerate(extended_input):
            if line.startswith("@SET") or line.startswith("@set"):
                components = line.split(None, 2)
                name = components[1]
                value = components[2]
                variables[name] = value
                logger.debug("Variable '{}' found with value '{}'".format(name, value))
            else:
                input_set_removed.append(line)

        # Place the variables
        variable_pattern = r"\$\{(\w+)\}|\$(\w+)"
        compiled = re.compile(variable_pattern)
        reserved = ("include", "set", "if", "endif")
        input_variables_replaced = []
        for line in input_set_removed:
            results = compiled.finditer(line)
            new_line = line
            offset = 0
            for result in results:
                options = result.groups()
                first = options[0]
                second = options[1]
                if first:
                    name = first
                elif second:
                    name = second
                if name in reserved:
                    continue
                value = variables.get(name)
                if not value:
                    logger.error("Value for variable '{}' not set.".format(name))
                    continue
                len_value = len(value)
                len_name = len(name)
                start = result.start()
                end = result.end()
                beginning = new_line[:offset+start]
                rest = new_line[offset+end:]
                new_line = beginning + value + rest
                offset += len_value - len_name - 1
            input_variables_replaced.append(new_line)

        self.input_lines = input_variables_replaced
