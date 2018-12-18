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

from builtins import next
from builtins import range
import os
import re
import logging
import importlib
from nomadcore.baseclasses import ParserInterface

# Needs to be imported in order for the importlib calls to work in python 2.7
import cp2kparser.versions.cp2k262.singlepointparser

logger = logging.getLogger("nomad")


class CP2KParser(ParserInterface):
    """This class handles the initial setup before any parsing can happen. It
    determines which version of CP2K was used to generate the output and then
    sets up a correct main parser.

    After the implementation has been setup, you can parse files with
    parse().
    """
    def __init__(self, metainfo_to_keep=None, backend=None, default_units=None, metainfo_units=None, debug=False, log_level=logging.ERROR, store=True):
        super(CP2KParser, self).__init__(metainfo_to_keep, backend, default_units, metainfo_units, debug, log_level, store)

    def setup_version(self):
        """Setups the version by looking at the output file and the version
        specified in it.
        """
        # Search for the CP2K version specification and the RUN_TYPE for the
        # calculation. The correct and optimized parser is initialized based on
        # this information.
        regex_version = re.compile(r" CP2K\| version string:\s+CP2K version ([\d\.]+)")
        regex_run_type = re.compile(r"\s+GLOBAL\| Run type\s+(.+)")
        n_lines = 100
        version_id = None
        run_type = None
        with open(self.parser_context.main_file, 'r') as outputfile:
            for i_line in range(n_lines):
                try:
                    line = next(outputfile)
                except StopIteration:
                    break
                result_version = regex_version.match(line)
                result_run_type = regex_run_type.match(line)
                if result_version:
                    version_id = result_version.group(1).replace('.', '')
                if result_run_type:
                    run_type = result_run_type.group(1)

        if version_id is None:
            msg = "Could not find a version specification from the given main file."
            logger.exception(msg)
            raise RuntimeError(msg)
        if run_type is None:
            msg = "Could not find a version specification from the given main file."
            logger.exception(msg)
            raise RuntimeError(msg)

        # Setup the root folder to the fileservice that is used to access files
        dirpath, filename = os.path.split(self.parser_context.main_file)
        dirpath = os.path.abspath(dirpath)
        self.parser_context.file_service.setup_root_folder(dirpath)
        self.parser_context.file_service.set_file_id(filename, "output")

        # Setup the correct main parser based on the version id. If no match
        # for the version is found, use the main parser for CP2K 2.6.2
        self.setup_main_parser({"version_id": version_id, "run_type": run_type})

    @staticmethod
    def get_mainfile_regex():
        regex_str = (
            "  \*\*\*\* \*\*\*\* \*\*\*\*\*\*  \*\*  PROGRAM STARTED AT\s.*\n"
            " \*\*\*\*\* \*\* \*\*\*  \*\*\* \*\*   PROGRAM STARTED ON\s*.*\n"
            " \*\*    \*\*\*\*   \*\*\*\*\*\*    PROGRAM STARTED BY .*\n"
            " \*\*\*\*\* \*\*    \*\* \*\* \*\*   PROGRAM PROCESS ID .*\n"
            "  \*\*\*\* \*\*  \*\*\*\*\*\*\*  \*\*  PROGRAM STARTED IN .*\n"
        )
        return regex_str

    def get_metainfo_filename(self):
        return "cp2k.nomadmetainfo.json"

    def get_parser_info(self):
        return {'name': 'cp2k-parser', 'version': '1.0'}

    def setup_main_parser(self, version_dictionary):
        """
        Setups a main parser class for this calculation. The main class can be
        different for each version and run type.

        Args:
            version_id: An integer representing the CP2K version. The version
                number is originally a string the form '2.6.2', but here the
                numbers are just concatenated into a single integer number 262.
            run_type: A string that identifies the RUN_TYPE for the
                calculation.  All the possible run types can be found in the
                CP2K reference manual.

        Returns:
            A python class that should be instantiated later with the correct
            parameters.
        """
        run_type = version_dictionary["run_type"]
        version_id = version_dictionary["version_id"]

        # Search for a RUN_TYPE specific parser
        parser_map = {
            "ENERGY": "SinglePointParser",
            "ENERGY_FORCE": "SinglePointParser",
            "WAVEFUNCTION_OPTIMIZATION": "SinglePointParser",
            "WFN_OPT": "SinglePointParser",
            "GEO_OPT": "GeoOptParser",
            "GEOMETRY_OPTIMIZATION": "GeoOptParser",
            "MD": "MDParser",
            "MOLECULAR_DYNAMICS": "MDParser",
        }
        try:
            parser = parser_map[run_type]
        except KeyError:
            logger.exception(
                "A parser corresponding to the run_type '{}' could not be found."
                .format(run_type)
            )
            raise

        # Currently the version id is a pure integer, so it can directly be mapped
        # into a package name.
        base = "cp2kparser.versions.cp2k{}.{}".format(version_id, parser.lower())
        parser_module = None
        parser_class = None

        try:
            parser_module = importlib.import_module(base)
        except ImportError:
            logger.warning(
                "Could not find a parser for version '{}' and run type '{}'. "
                "Trying to default to the base implementation for CP2K 2.6.2"
                .format(version_id, run_type)
            )
            base = "cp2kparser.versions.cp2k262.{}".format(parser.lower())
            try:
                parser_module = importlib.import_module(base)
            except ImportError:
                logger.exception(
                    "Tried to default to the CP2K 2.6.2 implementation but "
                    "could not find the correct modules for run_type '{}'."
                    .format(run_type)
                )
                raise
        try:
            parser_class = getattr(parser_module, "CP2K{}".format(parser))
        except AttributeError:
            logger.exception(
                "A parser class '{}' could not be found in the module '[]'."
                .format(parser_class, parser_module)
            )
            raise

        self.main_parser = parser_class(self.parser_context)
