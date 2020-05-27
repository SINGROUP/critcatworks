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

from builtins import str
from builtins import object
import numpy as np
import logging
from collections import defaultdict

logger = logging.getLogger("nomad")
metainfo_section_prefix = "x_cp2k_section_input_"
metainfo_data_prefix = "x_cp2k_input_"


class CP2KInput(object):
    """The contents of a CP2K simulation including default values and default
    units from the version-specific xml file.
    """

    def __init__(self, root_section):
        self.root_section = root_section

    @staticmethod
    def decode_cp2k_unit(unit):
        """Given a CP2K unit name, decode it as Pint unit definition.
        """
        map = {
            # Length
            "bohr": "bohr",
            "m": "meter",
            "pm": "picometer",
            "nm": "nanometer",
            "angstrom": "angstrom",

            # Angle
            "rad": "radian",
            "deg": "degree",

            #Energy
            "Ry": "rydberg"
        }
        pint_unit = map.get(unit)
        if pint_unit:
            return pint_unit
        else:
            logger.error("Unknown CP2K unit definition '{}'.".format(unit))

    def set_parameter(self, path, value):
        parameter, section = self.get_parameter_and_section(path)
        if section is None:
            message = "The CP2K input does not contain a section {}".format(path)
            logger.warning(message)
        if parameter is None:
            message = "The CP2K input section {} does not contain a SECTION_PARAMETER".format(path)
            logger.warning(message)
        else:
            parameter.value = value

    def set_keyword(self, path, value, full):
        keyword, section = self.get_keyword_and_section(path)
        # If keyword found, put data in there
        if keyword and section:
            keyword.value = value
        # Keyword not found in the input tree, assuming it is a default keyword
        elif section is not None:
            split_path = path.rsplit("/", 1)
            keyword = split_path[1]
            if section.default_keyword is not None:
                section.default_keyword.value += full + "\n"
            else:
                message = "The CP2K input does not contain the keyword {}, and there is no default keyword for the section {}".format(path, split_path[0])
                logger.warning(message)

    def get_section(self, path):
        split_path = path.split("/")
        section = self.root_section
        for part in split_path:
            section = section.get_subsection(part)
            if not section:
                message = "The CP2K input does not contain the section {}".format(path)
                logger.warning(message)
                return None
        return section

    def get_section_list(self, path):
        split_path = path.split("/")
        last_section = split_path[-1]
        split_path.pop()
        section = self.root_section
        for part in split_path:
            section = section.get_subsection(part)
            if not section:
                message = "The CP2K input does not contain the section {}".format(path)
                logger.warning(message)
                return None
        return section.get_subsection_list(last_section)

    def get_keyword_and_section(self, path):
        split_path = path.rsplit("/", 1)
        keyword = split_path[1]
        section_path = split_path[0]
        section = self.get_section(section_path)

        if section is None:
            message = "The CP2K input does not contain the section {}".format(path)
            logger.warning(message)
            return (None, None)

        keyword = section.get_keyword_object(keyword)
        if keyword and section:
            return (keyword, section)
        else:
            return (None, section)

    def get_keyword(self, path, raw=False, allow_default=True):

        keyword, section = self.get_keyword_and_section(path)
        if keyword:
            return keyword.get_value(raw, allow_default)

    def set_section_accessed(self, path):
        section = self.get_section(path)
        if section:
            section.accessed = True
        else:
            message = "The CP2K input does not contain the section {}".format(path)
            logger.warning(message)

    def get_default_unit(self, path):
        keyword, section = self.get_keyword_and_section(path)
        if keyword:
            return keyword.default_unit

    def get_unit(self, path):
        keyword, section = self.get_keyword_and_section(path)
        if keyword:
            return keyword.get_unit()

    def get_parameter_and_section(self, path):
        section = self.get_section(path)
        if section is None:
            return (None, None)
        if section.section_parameter is not None:
            parameter = section.section_parameter
            parameter.accessed = section.accessed
            return (parameter, section)
        else:
            return (None, section)

    def get_parameter(self, path):
        parameter, section = self.get_parameter_and_section(path)
        if parameter:
            if parameter.value:
                return parameter.value
            elif section and section.accessed:
                return parameter.lone_value


class Section(object):
    """An input section in a CP2K calculation.
    """
    __slots__ = ['accessed', 'name', 'keywords', 'default_keyword_names', 'default_keyword', 'section_parameter', 'sections', 'description']

    def __init__(self, name):
        self.accessed = False
        self.name = name
        self.keywords = defaultdict(list)
        self.default_keyword_names = []
        self.default_keyword = None
        self.section_parameter = None
        self.sections = defaultdict(list)
        self.description = None

    def get_keyword_object(self, name):
        keyword = self.keywords.get(name)
        if keyword:
            if len(keyword) == 1:
                return keyword[0]
            else:
                logger.error("The keyword '{}' in '{}' does not exist or has too many entries.".format(name, self.name))

    def get_keyword(self, name, raw=False, allow_default=True):
        """Returns the keyword value for the given name.

        Args:
            name: The name of the keyword
            raw: Boolean indicating if the raw value (not modified in any way)
                should be returned.
            allow_default: Boolean indicating if it is allowed to return the
                default value is no actual value was set by the user in the input.
        """
        keyword_object = self.get_keyword_object(name)
        return keyword_object.get_value(raw, allow_default)

    def get_subsection(self, name):
        subsection = self.sections.get(name)
        if subsection:
            if len(subsection) == 1:
                return subsection[0]
            else:
                logger.error("The subsection '{}' in '{}' has too many entries.".format(name, self.name))
        else:
            logger.error("The subsection '{}' in '{}' does not exist.".format(name, self.name))

    def get_subsection_list(self, name):
        subsection = self.sections.get(name)
        return subsection

    def get_section_parameter(self):
        """Get the section parameter, or if not specified the lone keyword
        value.
        """
        if self.section_parameter is not None:
            value = self.section_parameter.value
            if value is None:
                value = self.section_parameter.lone_keyword_value
            return value.upper()


class InputObject(object):
    """Base class for all kind of data elements in the CP2K input.
    """
    __slots__ = ['name', 'value', 'default_value', 'description', 'data_type', 'data_dimension']

    def __init__(self, name):
        self.name = name
        self.value = None
        self.description = None
        self.data_type = None
        self.data_dimension = None
        self.default_value = None


class Keyword(InputObject):
    """Information about a keyword in a CP2K calculation.
    """
    __slots__ = ['unit', 'value_no_unit', 'default_unit', 'default_name']

    def __init__(self, name, default_value,  default_unit, default_name):
        super(Keyword, self).__init__(name)
        self.unit = None
        self.value_no_unit = None
        self.default_unit = default_unit
        self.default_value = default_value
        self.default_name = default_name

    def get_value(self, raw=False, allow_default=True):
        if raw:
            return self._get_value_raw()
        else:
            return self._get_value_formatted(allow_default)

    def _get_value_raw(self):
        """Returns the unformatted value of this keyword. This is exactly what
        was set by the used in the input as a string.
        """
        return self.value

    def _get_value_formatted(self, allow_default=False):
        """Returns the value stored in this keyword by removing the possible
        unit definition and formatting the string into the correct data type.
        If asked, will use the default value if not actual value was set by
        user.
        """
        # Decode the unit and the value if not done before
        proper_value = None
        if self.default_unit is not None:
            if self.value_no_unit is None:
                self.decode_cp2k_unit_and_value()
        if self.value_no_unit is not None:
            proper_value = self.value_no_unit
        else:
            proper_value = self.value

        # if allow_default:
        if proper_value is None:
            proper_value = self.default_value
        if proper_value is None:
            return None

        returned = None
        dim = int(self.data_dimension)
        splitted = proper_value.split()
        if len(splitted) != dim:
            logger.error("The dimensions of the CP2K input parameter {} do not match the specification in the XML file.".format(self.name))

        if dim == 1:
            try:
                if self.data_type == "integer":
                    returned = int(proper_value)
                elif self.data_type == "real":
                    returned = float(proper_value)
                elif self.data_type == "word":
                    returned = str(proper_value)
                elif self.data_type == "keyword":
                    returned = str(proper_value)
                elif self.data_type == "string":
                    returned = str(proper_value)
                elif self.data_type == "logical":
                    returned = str(proper_value)
                else:
                    logger.error("Unknown data type '{}'".format(self.data_type))
                    return
            except TypeError:
                logger.error("The CP2K input parameter {} could not be converted to the type specified in the XML file.".format(self.name))
                return
        else:
            try:
                if self.data_type == "integer":
                    returned = np.array([int(x) for x in splitted])
                elif self.data_type == "real":
                    returned = np.array([float(x) for x in splitted])
                elif self.data_type == "word":
                    returned = np.array([str(x) for x in splitted])
                elif self.data_type == "keyword":
                    returned = np.array([str(x) for x in splitted])
                elif self.data_type == "string":
                    returned = np.array([str(x) for x in splitted])
                elif self.data_type == "logical":
                    returned = np.array([str(x) for x in splitted])
                else:
                    logger.error("Unknown data type '{}'".format(self.data_type))
                    return
            except TypeError:
                logger.error("The CP2K input parameter {} could not be converted to the type specified in the XML file.".format(self.name))
                return

        return returned

    def determine_value_and_unit(self):
        """If the units of this value can be changed, return a value and the
        unit separately.
        """
        if self.default_unit is not None:
            if self.value_no_unit is None:
                self.decode_cp2k_unit_and_value()
            return self.value_no_unit
        else:
            return self.value

    def get_unit(self):

        if self.unit is None:
            self.decode_cp2k_unit_and_value()
            if self.unit is not None:
                return self.unit
            elif self.default_unit is not None:
                return self.default_unit
        else:
            return self.unit

        logger.error("The keyword '{}' does not have a unit.".format(self.default_name))
        return None

    def decode_cp2k_unit_and_value(self):
        """Given a CP2K unit name, decode it as Pint unit definition.
        """
        if self.value is not None:
            splitted = self.value.split(None, 1)
            unit_definition = splitted[0]
            if unit_definition.startswith('[') and unit_definition.endswith(']'):
                unit_definition = unit_definition[1:-1]
                self.unit = CP2KInput.decode_cp2k_unit(unit_definition)
                self.value_no_unit = splitted[1]
            elif self.default_unit:
                logger.debug("No special unit definition found, returning default unit.")
                self.unit = CP2KInput.decode_cp2k_unit(self.default_unit)
                self.value_no_unit = self.value
            else:
                logger.debug("The value has no unit, returning bare value.")
                self.value_no_unit = self.value


class SectionParameters(InputObject):
    """Section parameters in a CP2K calculation.

    Section parameters are the short values that can be added right after a
    section name, e.g. &PRINT ON, where ON is the section parameter.
    """
    __slots__ = ['lone_keyword_value', 'accessed']

    def __init__(self, default_value, lone_keyword_value):
        super(SectionParameters, self).__init__("SECTION_PARAMETERS")
        self.default_value = default_value
        self.lone_keyword_value = lone_keyword_value
        self.accessed = None

    def get_value(self):
        """Returns the value for this section parameter. Uses the user given
        value primarily, if the section is not used, return the default value
        and if the section is defined but without explicit section parameter
        returns the lone keyword value.
        """
        if self.value is not None:
            return self.value
        elif self.accessed:
            return self.lone_keyword_value
        else:
            return self.default_value


class DefaultKeyword(InputObject):
    """Default keyword in the CP2K input.
    """
    __slots__ = ['lone_value']

    def __init__(self):
        super(DefaultKeyword, self).__init__("DEFAULT_KEYWORD")
        self.lone_value = None
        self.value = ""
