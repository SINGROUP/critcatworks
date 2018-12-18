from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from enum import Enum
from operator import itemgetter
from builtins import range
from builtins import object
from builtins import str
import numpy as np
import logging
from collections import defaultdict
logger = logging.getLogger(__name__)


class ParserEvent(Enum):
    """Enumerations for the different parser events when traversing the
    results.
    """
    open_section = 1
    close_section = 2
    add_value = 3
    add_array_value = 4


class ParserKeyError(Exception):
    """A custom exception for all cases where a certain metainfo can not be
    found in the results.
    """
    pass


class DummyFile(object):
    """Mimics a file object by defininf a write interface, but actually does
    not write anything. This allows one to used this object in functions which
    require a file object as input, but you dont want to actually write
    anything.
    """
    def write(self, input):
        pass


class LocalBackend(object):
    """A backend that outputs results into a regular python dictionary. This is
    useful if you wan't to run the parser with python only.
    """
    def __init__(self, metaInfoEnv, debug=True, store=True):
        """
        Args:
            metaInfoEnv: The description of the metainfo environment.
            debug: Boolean indicating whether some debugging should be done
                (check that correct types and shapes are pushed into backend)
            store: Boolean indicating whether the parsed results should be
                stored in memory. Useful to skip the storing if you just want
                to quickly check that everything is runnning fine and don't
                want to waste RAM in doing so.
        """
        self.__metaInfoEnv = metaInfoEnv
        self.fileOut = DummyFile()
        self.__gIndex = -1
        self.__openSections = set()
        self.__lastIndex = {}
        self.stats = {}
        self.dataManagers = {}
        self.sectionManagers = {}
        self.results = Results(metaInfoEnv, self.dataManagers, self.sectionManagers)
        self.debug = debug
        self.store = store

        for ikNames, ik in metaInfoEnv.infoKinds.items():
            if ik.kindStr == "type_section":
                parentS = list(metaInfoEnv.firstAncestorsByType(ik.name).get("type_section", [[]])[0])
                parentS.sort()
                self.sectionManagers[ik.name] = SectionManager(
                    metaInfo=ik,
                    parentSectionNames=parentS, debug=self.debug)
        for ikNames, ik in metaInfoEnv.infoKinds.items():
            if ik.kindStr == "type_document_content" or ik.kindStr == "type_dimension":
                superSectionNames = metaInfoEnv.firstAncestorsByType(ik.name).get("type_section", [[]])[0]
                if not superSectionNames:
                    raise Exception("MetaInfo of conrete value %s is not in any superSection" % ik.name)
                elif len(superSectionNames) > 1:
                    raise Exception("MetaInfo of concrete value %s has multiple superSections (%s)" %
                                    (ik.name, superSectionNames))
                self.dataManagers[ik.name] = DataManager(ik, self.sectionManagers[superSectionNames[0]])

    def openSection(self, metaName):
        """opens a new section and returns its new unique gIndex"""
        manager = self.sectionManagers[metaName]
        newIndex = manager.openSection(self)
        self.__openSections.add((metaName, newIndex))
        return newIndex

    def openNonOverlappingSection(self, metaName):
        """opens a new non overlapping section"""
        if any(x[0] == metaName for x in self.__openSections):
            raise Exception("Section %s is not supposed to overlap" % metaName)
        return self.openSection(metaName)

    def closeSection(self, metaName, gIndex):
        manager = self.sectionManagers[metaName]
        manager.closeSection(self, gIndex)
        if (metaName, gIndex) in self.__openSections:
            self.__openSections.remove((metaName, gIndex))

    def closeNonOverlappingSection(self, metaName):
        """closes a non overlapping section"""
        openGIndexes = [x for x in self.__openSections if x[0] == metaName]
        if len(openGIndexes) != 1:
            if not openGIndexes:
                raise Exception("Call to closeNonOverlapping(%s) with no open section" % metaName)
            else:
                raise Exception("Section %s was not supposed to overlap, found %s open when closing" % (metaName, openGIndexes))
        self.closeSection(metaName, openGIndexes[0][1])

    def addValue(self, metaName, value, gIndex=-1):

        dataManager = self.dataManagers[metaName]

        if self.debug:

            # Check that the value is actually of scalar type
            value_type = type(value)
            if value_type not in [float, int, bool, type(b""), type(u""), str, np.float64]:
                raise TypeError("Could not use function 'addValue' to push value '{}' with type '{}' for metainfo '{}'.".format(value, value_type, metaName))

            # Check that the metainfo should be scalar
            metainfo_shape = dataManager.metaInfo.shape
            if metainfo_shape is not None:
                if len(metainfo_shape) != 0:
                    raise TypeError("The metainfo '{}' does not support scalar values. Check the shape attribute of the metainfo and use the function addArrayValues() instead if the value should be an array.".format(metaName))

            # Check the type
            dtype_str = dataManager.metaInfo.dtypeStr
            single_types = self.single_value_type_for_metainfo_type(dtype_str)
            actual_numpy_type = type(value)
            if actual_numpy_type not in single_types:
                raise TypeError("The given value for metainfo '{}' is of incorrrect type. The type was '{}' when it should be one of '{}'".format(metaName, actual_numpy_type, single_types))

        dataManager.superSectionManager.addValue(dataManager.metaInfo, value, gIndex)

    def addRealValue(self, metaName, value, gIndex=-1):
        self.addValue(metaName, value, gIndex)

    def addArrayValues(self, metaName, values, gIndex=-1):

        dataManager = self.dataManagers[metaName]

        if self.debug:

            # Check that the value is actually a numpy array
            if not isinstance(values, np.ndarray):
                raise TypeError("The value provided for '{}' is not a valid numpy array. Please only push numpy arrays with the backend function addArrayValues().".format(metaName))

            # Check that the metainfo should be an array
            metainfo_shape = dataManager.metaInfo.shape
            if len(metainfo_shape) == 0:
                raise TypeError("The metainfo '{}' does not support arrays. Check the shape attribute of the metainfo and use the function addValue() instead if the value should be scalar.".format(metaName))

            # Check the number of dimensions
            array_shape = values.shape
            len_meta_dim = len(metainfo_shape)
            len_array_dim = len(array_shape)
            if len_array_dim != len_meta_dim:
                raise TypeError("Incompatible shape provided for metainfo '{}'. The shape was '{}' whereas it should be '{}'. Check the shape attribute of the metainfo".format(metaName, array_shape, metainfo_shape))

            # If the shapes are given as integers in the metainfo we can also
            # check the number of values in each dimension
            try:
                [int(x) for x in metainfo_shape]
            except Exception:
                pass
            else:
                for index in range(len_meta_dim):
                    array_dim = array_shape[index]
                    metainfo_dim = metainfo_shape[index]
                    if array_dim != metainfo_dim:
                        raise TypeError("Incompatible shape provided for metainfo '{}'. The shape was '{}' whereas it should be '{}'. Check the shape attribute of the metainfo".format(metaName, array_shape, metainfo_shape))

            # Check the type
            dtype_str = dataManager.metaInfo.dtypeStr
            array_types = self.array_type_for_metainfo_type(dtype_str)
            actual_numpy_type = values.dtype.type
            if actual_numpy_type not in array_types:
                raise TypeError("The given array for metainfo '{}' has incorrect type of values in it. The values given are '{}', whereas the datatype given in metainfo is '{}'".format(metaName, actual_numpy_type, dtype_str))

        dataManager.superSectionManager.addArrayValues(dataManager.metaInfo, values, gIndex)

    def array_type_for_metainfo_type(self, dtypeStr):
        """Returns a list of numpy types correspoding to the dtypeStr of a
        metainfo.
        """
        if dtypeStr == "f":
            return [np.float_, np.float64, np.float32, np.int_, np.int64, np.int32, np.int16, np.int8]
        elif dtypeStr == "i":
            return [np.int_, np.int64, np.int32, np.int16, np.int8]
        elif dtypeStr == "b":
            return [np.bool_]
        elif dtypeStr == "C":
            return [np.string_, np.unicode_]
        elif dtypeStr == "r":
            return [np.int_, np.int64, np.int32, np.int16, np.int8]
        else:
            raise TypeError("Could not determine the numpy type for metainfo type '{}'".format(dtypeStr))

    def single_value_type_for_metainfo_type(self, dtypeStr):
        """Returns a list of numpy types corresponding to the dtypeStr of a
        metainfo.
        """
        if dtypeStr == "f":
            return [float, int, np.float_, np.float64, np.float32]
        elif dtypeStr == "i":
            return [int, np.int_, np.int64, np.int32, np.int16, np.int8]
        elif dtypeStr == "b":
            return [bool, np.bool_]
        elif dtypeStr == "C":
            return [type(b""), type(u""), str, np.string_, np.unicode_]
        elif dtypeStr == "r":
            return [int, np.int_, np.int64, np.int32, np.int16, np.int8]
        else:
            raise TypeError("Could not determine the type for metainfo type '{}'".format(dtypeStr))

    def setArrayValues(self, metaName, values, offset=None, gIndex=-1, unit=None):
        """Adds values to the last array added, array must be a numpy array
        """
        dataManager = self.dataManagers[metaName]
        dataManager.superSectionManager.setArrayValues(dataManager.metaInfo, values, offset, gIndex)

    def metaInfoEnv(self):
        return self.__metaInfoEnv

    def startedParsingSession(self, mainFileUri, parserInfo, parserStatus=None, parserErrors=None):
        pass

    def finishedParsingSession(self, parserStatus, parserErrors, mainFileUri=None, parserInfo=None,
                               parsingStats=None):
        """Called when the parsing finishes.
        """
        pass

    def addMatchTelemetry(self, match_telemetry, gIndex=-1):
        """ should be called for outputting match telemetry data:
        input data, together with capture info """
        pass

    def pwarn(self, msg):
        """Used to catch parser warnings. Currently disabled in the local
        backend.
        """
        pass


#===============================================================================
class Results(object):
    """A wrapper object for the collection of results gathered by a parser.
    """
    def __init__(self, metaInfoEnv, datamanagers, sectionmanagers):
        self._datamanagers = datamanagers
        self._sectionmanagers = sectionmanagers
        self._shortnames = defaultdict(list)
        self._metaInfoEnv = metaInfoEnv

    def __getitem__(self, metaname):
        """Return the data or section corrresponding the the given metainfo
        name. If given a section name, this function will return a list of
        Section objects. If given a name of a concrete value, this function
        will return all instances of that value as a list.

        Args:
            metaname: The unique name of the metainfo to get.

        Raises:
            LookupError: if the metaname is not defined in the metainfo
                environment or the parser has not output any value for it.
            ParserKeyError: if the parser did not output the queried metainfo.
        """
        self.test_validity(metaname)

        # See if in sections
        sectionmanager = self._sectionmanagers.get(metaname)
        if sectionmanager is not None:
            return sectionmanager.openSections

        # See if in data
        datamanager = self._datamanagers.get(metaname)
        if datamanager is not None:
            sectionmanager = datamanager.superSectionManager
            open_sections = sectionmanager.openSections
            result = []
            for section in open_sections:
                try:
                    data = section[metaname]
                except KeyError:
                    pass
                else:
                    result.append(data)
            if len(result) == 1:
                return result[0]
            elif len(result) == 0:
                raise KeyError("Could not find a parsing result for '{}'. The parser did not output this value.".format(metaname))
            else:
                return np.array(result)

        raise LookupError("The metainfo definition doesn't seem to contain '{}'. Check for typos of update you metainfo repository.".format(metaname))

    def test_validity(self, metaname):
        """Tests if the given metaname is present in the metainfo environment.
        """
        metainfo = self._metaInfoEnv.infoKinds.get(metaname)
        if metainfo is None:
            raise LookupError("The metainfo name '{}' does not exist in the metainfo environment. Check for typos or try updating the metainfo git package.".format(metaname))

    def traverse(self):
        """A generator function for traversing the data in the parser results.

        This generator returns a tuple of three item: the metainfo name, the
        event type, and the event value.
        """
        root = self._sectionmanagers["section_run"]
        for x in self.traverse_recursive("section_run", root.openSections):
            yield x

    def traverse_recursive(self, name, open_sections):
        """A generator function for traversing the data in the parser results.
        """
        for i_section, section in enumerate(open_sections):
            yield (name, ParserEvent.open_section, i_section)

            key_to_type_map = {}
            simple_keys = list(section.simple_values.keys())
            for key in simple_keys:
                key_to_type_map[key] = "simple"
            array_keys = list(section.array_values.keys())
            for key in array_keys:
                key_to_type_map[key] = "array"
            subsection_keys = list(section.subsections.keys())
            for key in subsection_keys:
                key_to_type_map[key] = "subsection"

            keys = []
            keys.extend(simple_keys)
            keys.extend(array_keys)
            keys.extend(subsection_keys)
            keys = sorted(keys)

            for key in keys:
                key_type = key_to_type_map[key]
                if key_type == "simple":
                    simple_value = section.simple_values[key]
                    yield (key, ParserEvent.add_value, simple_value)
                elif key_type == "array":
                    array_value = section.array_values[key]
                    yield (key, ParserEvent.add_array_value, array_value)
                elif key_type == "subsection":
                    subsection_value = section.subsections[key]
                    for x in self.traverse_recursive(key, subsection_value):
                        yield x
                else:
                    raise KeyError("Trying to access unknown data type.")

            yield (name, ParserEvent.close_section, i_section)

            # for value_name, value_value in section.simple_values.items():
                # yield (value_name, ParserEvent.add_value, value_value)
            # for array_name, array_value in section.array_values.items():
                # yield (array_name, ParserEvent.add_array_value, array_value)
            # for x in self.traverse_recursive(section.subsections):
                # yield x

    def print_summary(self):
        """Return a string representing the data contained in the results. This
        is a summary that can be used for debugging.
        """
        metas = {}
        roots = {}

        for meta in self._metaInfoEnv.infoKinds.values():
            metaobj = {}
            metaobj["name"] = meta.name
            metaobj["children"] = []
            metaobj["parents"] = meta.superNames
            metaobj["kindStr"] = meta.kindStr
            mapping = {
                "type_section": 0,
                "type_abstract_document_content": 1,
                "type_document_content": 2,
                "type_dimension": 3,
                "type_meta": 4,
            }
            metaobj["kind_number"] = mapping.get(meta.kindStr)
            metas[meta.name] = metaobj

        for meta in metas.values():
            parentNames = meta["parents"]
            if len(parentNames) == 0:
                roots[meta["name"]] = meta
            else:
                for parentName in parentNames:
                    parent = metas[parentName]
                    parent["children"].append(meta)

        # Sort the children according to type
        for meta in metas.values():
            meta["children"].sort(key=itemgetter('kind_number', 'name'))

        section_run = roots["section_run"]
        self.print_metainfo(section_run)

    def print_metainfo(self, meta, level=0):
        """Recursive printing function for the metainfos. To print the whole
        tree, call this function on the root section.
        """
        name = meta["name"]
        metatype = meta["kindStr"]
        if metatype != "type_abstract_document_content":
            try:
                result = self[name]
            except LookupError:
                return
            if isinstance(result, dict):
                if len(result.keys()) == 0:
                    return

            if metatype == "type_section":
                print(level*"  " + name + ":")
            elif metatype == "type_document_content":
                print(level*"  " + name)
            elif metatype == "type_dimension":
                print(level*"  " + name)
            level += 1

        for child in meta["children"]:
            self.print_metainfo(child, level)


class Section(object):
    """Represents an open section.
    """
    def __init__(self, gIndex, references, parents, name, backend, debug=True):
        self.gIndex = gIndex
        self.references = references
        self.simple_values = {}
        self.array_values = {}
        self.subsections = {}
        self.parents = parents
        self.name = name
        self.backend = backend
        self.debug = debug
        self.has_results = False

    def __getitem__(self, metaName):
        """Returns the cached values corresponding to metaName. You can search
        values and subsections.
        """
        res = self.simple_values.get(metaName, None)
        if res is not None:
            return res
        res = self.array_values.get(metaName, None)
        if res is not None:
            return res
        res = self.subsections.get(metaName, None)
        if res is not None:
            return res

        raise KeyError(
            "The metainfo '{}' could not be found in the section '{}' with gIndex '{}'"
            .format(metaName, self.name, self.gIndex))

    def get(self, key, default=None):
        try:
            value = self[key]
        except KeyError:
            return default
        else:
            return value

    def __len__(self):
        n_simple_values = len(self.simple_values)
        n_array_values = len(self.array_values)
        n_subsections = len(self.subsections)
        return n_simple_values + n_array_values + n_subsections

    def keys(self):
        keys_simple_values = self.simple_values.keys()
        keys_array_values = self.array_values.keys()
        keys_subsections = self.subsections.keys()

        for key in keys_simple_values:
            yield key
        for key in keys_array_values:
            yield key
        for key in keys_subsections:
            yield key

    def values(self):
        values_simple_values = self.simple_values.values()
        values_array_values = self.array_values.values()
        values_subsections = self.subsections.values()

        for value in values_simple_values:
            yield value
        for value in values_array_values:
            yield value
        for value in values_subsections:
            yield value

    def items(self):
        items_simple_values = self.simple_values.items()
        items_array_values = self.array_values.items()
        items_subsections = self.subsections.items()

        for item in items_simple_values:
            yield item
        for item in items_array_values:
            yield item
        for item in items_subsections:
            yield item

    def __contains__(self, key):
        keys = self.keys()
        return key in keys

    def addValue(self, metaInfo, value):
        if self.backend.store:
            if self.debug:
                vals = self.simple_values.get(metaInfo.name, None)
                if vals is None:
                    self.simple_values[metaInfo.name] = value
                else:
                    raise Exception("Trying to add values multiple times for metaname {} in section {}. ".format(metaInfo.name, self.name))
            else:
                self.simple_values[metaInfo.name] = value

    def setArrayValues(self, metaInfo, values, offset=None):
        if self.backend.store:
            vals = self.array_values.get(metaInfo.name, None)
            if vals is None:
                raise Exception("setArrayValues(%s,...) called before adding a value" % metaInfo.name)
            else:
                if offset:
                    idxs = [slice(offset[i], offset[i] + values.shape[i]) for i in range(len(offset))]
                else:
                    idxs = [slice(0, x) for x in values.shape]
                vals[len(vals) - 1][idxs] = values

    def addArrayValues(self, metaInfo, values):
        if self.backend.store:
            if self.debug:
                vals = self.array_values.get(metaInfo.name, None)
                if vals is None:
                    self.array_values[metaInfo.name] = values
                else:
                    raise Exception("Trying to add values multiple times for metaname {} in section {}. ".format(metaInfo.name, self.name))
            else:
                self.array_values[metaInfo.name] = values

    def addSubsection(self, metaInfo, section):
        vals = self.subsections.get(metaInfo.name, None)
        if vals is None:
            self.subsections[metaInfo.name] = [section]
        else:
            vals.append(section)


class SectionManager(object):
    """Manages the sections for the given metainfo.
    """
    def __init__(self, metaInfo, parentSectionNames, lastSectionGIndex=-1, debug=True):
        self.metaInfo = metaInfo
        self.parentSectionNames = parentSectionNames
        self.lastSectionGIndex = lastSectionGIndex
        self.debug = debug
        self.openSections = []

    def openSection(self, backend):
        newGIndex = self.lastSectionGIndex + 1
        self.openSectionWithGIndex(backend, newGIndex)
        return newGIndex

    def openSectionWithGIndex(self, backend, gIndex):
        self.lastSectionGIndex = gIndex
        references = []
        parents = []
        parent_found = False
        for parentName in self.parentSectionNames:
            pSect = backend.sectionManagers.get(parentName)
            try:
                parentSection = pSect.openSections[pSect.lastSectionGIndex]
            except KeyError:
                pass
            else:
                parent_found = True
                parents.append(parentSection)
            if pSect:
                references.append(pSect.lastSectionGIndex)
            else:
                references.append(-1)

        # If the section is supposed to have parents, and none were actually
        # open, raise an error
        if not parent_found and len(self.parentSectionNames) != 0:
            raise LookupError("Could not open section '{}' because none of it's parent sections '{}' could not be found".format(self.metaInfo.name, self.parentSectionNames))

        new_section = Section(gIndex, references, parents, self.metaInfo.name, backend, debug=self.debug)
        self.openSections.append(new_section)
        if parent_found:
            parents[0].addSubsection(self.metaInfo, new_section)

    def closeSection(self, backend, gIndex):
        pass

    def addValue(self, valueMetaInfo, value, gIndex):
        if (gIndex == -1):
            gI = self.lastSectionGIndex
        else:
            gI = gIndex
        try:
            self.openSections[gI].addValue(valueMetaInfo, value)
        except KeyError:
            raise Exception("Cannot add value for metadata %s to section %d (%d) of %s, as it is not open" % (valueMetaInfo.name, gI, gIndex, self.metaInfo.name))

    def setArrayValues(self, valueMetaInfo, value, offset=None, gIndex=-1):
        if gIndex == -1:
            gI = self.lastSectionGIndex
        else:
            gI = gIndex
        try:
            self.openSections[gI].setArrayValues(valueMetaInfo, value, offset)
        except KeyError:
            raise Exception("Cannot set array values for metadata %s to section %d (%d) of %s, as it is not open" % (valueMetaInfo.name, gI, gIndex, self.metaInfo.name))

    def addArrayValues(self, valueMetaInfo, value, offset=None, gIndex=-1):
        if gIndex == -1:
            gI = self.lastSectionGIndex
        else:
            gI = gIndex
        try:
            self.openSections[gI].addArrayValues(valueMetaInfo, value)
        except KeyError:
            raise Exception("Cannot add array values for metadata %s to section %d (%d) of %s, as it is not open" % (valueMetaInfo.name, gI, gIndex, self.metaInfo.name))


class DataManager(object):
    """Stores the parent (SectionManager) for the given metainfo.
    """
    def __init__(self, metaInfo, superSectionManager):
        self.metaInfo = metaInfo
        self.superSectionManager = superSectionManager
