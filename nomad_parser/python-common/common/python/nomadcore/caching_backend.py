from builtins import str
from builtins import zip
from builtins import range
from builtins import object
from functools import reduce
import json
import numpy as np
import sys
import logging
from nomadcore import parser_backend
from nomadcore.unit_conversion import unit_conversion

class CachingLevel(object):
    Forward = 1
    Cache = 2
    ForwardAndCache = 3
    Ignore = 4
    PreOpenedCache = 5
    PreOpenedIgnore = 6

    @classmethod
    def toText(cls, cl):
        if cl == cls.Forward:
            return "Forward"
        elif cl == cls.ForwardAndCache:
            return "ForwardAndCache"
        elif cl == cls.Cache:
            return "Cache"
        elif cl == cls.Ignore:
            return "Ignore"
        elif cl == cls.PreOpenedIgnore:
            return "PreOpenedIgnore"
        elif cl == cls.PreOpenedCache:
            return "PreOpenedIgnore"
        else:
            return "CachingLevel(%s)" % cl

    @classmethod
    def restrict(cls, a, b):
        """restricts a keeping in account that it is contained in b (i.e. if b is not forwarded then also a cannot be forwarded)"""
        if b == cls.Cache or b == cls.Ignore:
            if a == cls.Forward:
                return cls.Ignore
            elif a == cls.ForwardAndCache:
                return cls.Cache
        return a
class CachingSection(object):
    """Represents an open section, and can cache values"""
    def __init__(self, gIndex, references, storeInSuper = False, activeBackend=None):
        self.gIndex = gIndex
        self.references = references
        self.simpleValues = {}
        self.arrayValues = {}
        self.subSectionValues = {}
        self.storeInSuper = storeInSuper
        self.activeBackend = activeBackend

    def __getitem__(self, metaName):
        """Returns the cached values corresponding to metaName, or None if
        there aren't any. You can search values and subsections.
        """
        res = self.simpleValues.get(metaName, None)
        if res:
            return res
        res = self.arrayValues.get(metaName, None)
        if res:
            return res
        res = self.subSectionValues.get(metaName, None)
        return res

    def get_latest_value(self, metaname_list):
        """Goes through the given list of metanames and tries to find the data
        by succesfully opening the latest value from within this section and
        it's subsections.  If the given path leads to an actual value, it is
        returned.

        Args:
            metaname_list: list of metainfo names

        Returns:
            The value that if found when successively opening the data in this
            section and it's subsection. If some of the metanames in the given
            list are not found, returns None.
        """
        # Make path into a list
        if not isinstance(metaname_list, (tuple, list)):
            metaname_list = [metaname_list]

        # Go deeper into the section structure by iterating through the given
        # list of metanames.
        current = self
        for metaname in metaname_list:
            new_section = current[metaname]
            if new_section is not None:
                current = new_section[-1]
            else:
                return None

        return current

    def add_latest_value(self, metaname_list, metaname, require=False):
        latest_value = self.get_latest_value(metaname_list)
        if latest_value is not None:
            self.activeBackend.addValue(metaname, latest_value)
        else:
            if require:
                raise KeyError("Could not add the latest value from path '{}' because not value was found.")

    def add_latest_array_values(self, metaname_list, metaname, require=False):
        latest_value = self.get_latest_value(metaname_list)
        if latest_value is not None:
            self.activeBackend.addArrayValues(metaname, latest_value)
        else:
            if require:
                raise KeyError("Could not add the latest value from path '{}' because not value was found.")

    def addValue(self, metaInfo, value):
        vals = self.simpleValues.get(metaInfo.name, None)
        if vals is None:
            self.simpleValues[metaInfo.name] = [value]
        else:
            vals.append(value)

    def setArrayValues(self, metaInfo, values, offset = None):
        vals = self.arrayValues.get(metaInfo.name, None)
        if vals is None:
            raise Exception("setArrayValues(%s,...) called before adding a value" % metaInfo.name)
        else:
            if offset:
                idxs = [slice(offset[i], offset[i] + values.shape[i]) for i in range(len(offset))]
            else:
                idxs = [slice(0, x) for x in values.shape]
            vals[len(vals) - 1][idxs] = values

    def addArrayValues(self, metaInfo, values):
        vals = self.arrayValues.get(metaInfo.name, None)
        if vals is None:
            self.arrayValues[metaInfo.name] = [values]
        else:
            vals.append(values)

    def addSubsection(self, metaInfo, section):
        vals = self.subSectionValues.get(metaInfo.name, None)
        if vals is None:
            self.subSectionValues[metaInfo.name] = [section]
        else:
            vals.append(section)

    def __str__(self):
        return """CachingSection{{
        gIndex:{0},
        references:{1},
        simpleValues:{2},
        arrayValues.keys:{3},
        subSectionValues.keys:{4},
        storeInSuper:{5},
    }}""".format(
        self.gIndex,
        self.references,
        self.simpleValues,
        list(self.arrayValues.keys()),
        list(self.subSectionValues.keys()),
        self.storeInSuper)

class CachingSectionManager(object):

    def __init__(self, metaInfo, parentSectionNames, lastSectionGIndex = -1, openSections = {}, onClose = [], onOpen = [], storeInSuper = False, forwardOpenClose = True, preOpened = False):
        self.metaInfo = metaInfo
        self.parentSectionNames = parentSectionNames
        self.lastSectionGIndex = lastSectionGIndex
        self.openSections = dict(openSections)
        self.onClose = list(onClose)
        self.onOpen = list(onOpen)
        self.storeInSuper = storeInSuper
        self.forwardOpenClose = forwardOpenClose
        self.preOpened = preOpened

    def setSectionInfo(self, gIndex, references):
        self.openSections[gIndex].references = [references[x] for x in self.parentSectionNames]

    def get_latest_section(self):
        """Returns the latest opened section if it is still open.
        """
        section = self.openSections.get(self.lastSectionGIndex)
        return section

    def openSection(self, backend):
        newGIndex = self.lastSectionGIndex + 1
        self.openSectionWithGIndex(backend, newGIndex)
        return newGIndex

    def openSectionWithGIndex(self, backend, gIndex):
        self.lastSectionGIndex = gIndex
        references = []
        for parentName in self.parentSectionNames:
            pSect = backend.sectionManagers.get(parentName, None)
            if pSect:
                references.append(pSect.lastSectionGIndex)
            else:
                references.append(-1)
        opened = CachingSection(gIndex, references, storeInSuper = self.storeInSuper, activeBackend=backend)
        self.openSections[gIndex] = opened
        for onOpenF in self.onOpen:
            onOpenF(backend, gIndex, opened)

    def closeSection(self, backend, gIndex):
        toClose = self.openSections[gIndex]
        for onCloseF in self.onClose:
            onCloseF(backend, gIndex, toClose)
        if toClose.storeInSuper:
            for superGIndex, superName in zip(toClose.references,self.parentSectionNames):
                superSect = backend.sectionManagers[superName].openSections.get(superGIndex, None)
                if superSect:
                    superSect.addSubsection(self.metaInfo, toClose)
                else:
                    backend.storeToClosedSuper(superName, superGIndex, self.metaInfo, toClose)
            if not toClose.references: # checks for empty list
                backend.parsingSession.addSubsection(self.metaInfo, toClose)
        if self.forwardOpenClose and backend.superBackend:
            backend.superBackend.closeSection(self.metaInfo.name, gIndex)
        del self.openSections[gIndex]

    def openSectionInfo(self, gIndex):
        section = self.openSections.get(gIndex, None)
        if section:
            return "section {name} ({references})".format(
                name = section.name,
                references = list(zip(self.parentSectionNames, section.references)))
        else:
            return "section {gIndex} in {name} is not open!!".format(
                name = self.metaInfo.name,
                gIndex = gIndex)

    def addValue(self, valueMetaInfo, value, gIndex):
        if (gIndex == -1):
            gI = self.lastSectionGIndex
        else:
            gI = gIndex
        try:
            self.openSections[gI].addValue(valueMetaInfo, value)
        except:
            raise Exception("Cannot add value for metadata %s to section %d (%d) of %s, as it is not open" % (valueMetaInfo.name, gI, gIndex, self.metaInfo.name))

    def setArrayValues(self, valueMetaInfo, value, offset = None, gIndex = -1):
        if gIndex == -1:
            gI = self.lastSectionGIndex
        else:
            gI = gIndex
        try:
            self.openSections[gI].setArrayValues(valueMetaInfo, value, offset)
        except:
            raise Exception("Cannot set array values for metadata %s to section %d (%d) of %s, as it is not open" % (valueMetaInfo.name, gI, gIndex, self.metaInfo.name))

    def addArrayValues(self, valueMetaInfo, value, offset = None, gIndex = -1):
        if gIndex == -1:
            gI = self.lastSectionGIndex
        else:
            gI = gIndex
        try:
            self.openSections[gI].addArrayValues(valueMetaInfo, value)
        except:
            raise Exception("Cannot add array values for metadata %s to section %d (%d) of %s, as it is not open" % (valueMetaInfo.name, gI, gIndex, self.metaInfo.name))

class CachingDataManager(object):
    def __init__(self, metaInfo, superSectionManager, cachingLevel):
        self.metaInfo = metaInfo
        self.superSectionManager = superSectionManager
        self.cachingLevel = cachingLevel

class ActiveBackend(object):
    def __init__(self, metaInfoEnv, sectionManagers, dataManagers, superBackend, propagateStartFinishParsing = True, default_units=None, metainfo_units=None):
        self.__metaInfoEnv = metaInfoEnv
        self.sectionManagers = sectionManagers
        self.dataManagers = dataManagers
        self.superBackend = superBackend
        self.propagateStartFinishParsing = propagateStartFinishParsing
        self.default_units = default_units  # A mapping between dimension and an unit definition.
        self.metainfo_units = metainfo_units  # A mapping between metaname and an unit definition.

    @classmethod
    def activeBackend(cls, metaInfoEnv, cachingLevelForMetaName = {}, defaultDataCachingLevel = CachingLevel.ForwardAndCache, defaultSectionCachingLevel = CachingLevel.Forward, superBackend = None,
                      onClose = {}, onOpen = {}, propagateStartFinishParsing = True, default_units=None, metainfo_units=None):
        for sectionName in onClose.keys():
            if not sectionName in metaInfoEnv:
                raise Exception("Found trigger for non existing section %s" % sectionName)
            elif metaInfoEnv.infoKinds[sectionName].kindStr != "type_section":
                raise Exception("Found trigger for %s which is not a section but %s" %
                                (sectionName, json.dumps(metaInfoEnv.infoKinds[sectionName].toDict(), indent=2)))
        for sectionName in onOpen.keys():
            if not sectionName in metaInfoEnv:
                raise Exception("Found trigger for non existing section %s" % sectionName)
            elif metaInfoEnv.infoKinds[sectionName].kindStr != "type_section":
                raise Exception("Found trigger for %s which is not a section but %s" %
                                (sectionName, json.dumps(metaInfoEnv.infoKinds[sectionName].toDict(), indent=2)))
        sectionManagers = {}
        for ikNames, ik in metaInfoEnv.infoKinds.items():
            if ik.kindStr == "type_section":
                parentS, parentO = list(metaInfoEnv.firstAncestorsByType(ik.name).get("type_section", [[],[]]))
                parentS.sort()
                cachingLevel = reduce(CachingLevel.restrict, [cachingLevelForMetaName.get(x, defaultSectionCachingLevel) for x in ([ik.name] + parentS + parentO)])
                sectionManagers[ik.name] = CachingSectionManager(
                    metaInfo           = ik,
                    parentSectionNames = parentS,
                    storeInSuper       = (cachingLevel == CachingLevel.ForwardAndCache or cachingLevel == CachingLevel.Cache or cachingLevel == CachingLevel.PreOpenedCache),
                    forwardOpenClose   = (cachingLevel == CachingLevel.Forward or cachingLevel == CachingLevel.ForwardAndCache),
                    preOpened = (cachingLevel == CachingLevel.PreOpenedCache or cachingLevel == CachingLevel.PreOpenedIgnore),
                    onClose = onClose.get(ik.name, []),
                    onOpen = onOpen.get(ik.name, []))
        dataManagers = {}
        for ikNames, ik in metaInfoEnv.infoKinds.items():
            if ik.kindStr == "type_document_content" or ik.kindStr == "type_dimension":
                superSectionNames = metaInfoEnv.firstAncestorsByType(ik.name).get("type_section", [[]])[0]
                if not superSectionNames:
                    raise Exception("MetaInfo of conrete value %s is not in any superSection" % ik.name)
                elif len(superSectionNames) > 1:
                    raise Exception("MetaInfo of concrete value %s has multiple superSections (%s)" %
                                    (ik.name, superSectionNames))
                sectionManager = sectionManagers[superSectionNames[0]]
                dataManagers[ik.name] = CachingDataManager(ik, sectionManager,
                                                           CachingLevel.restrict(cachingLevelForMetaName.get(ik.name, defaultDataCachingLevel), CachingLevel.Forward if sectionManager.forwardOpenClose or sectionManager.preOpened else CachingLevel.Ignore))
        return ActiveBackend(metaInfoEnv, sectionManagers, dataManagers, superBackend, propagateStartFinishParsing, default_units, metainfo_units)

    def appendOnClose(self, sectionName, onClose):
        self.sectionManagers.onClose.append(onClose)

    def get_latest_section(self, metaname):
        """Returns the latest opened section with the given metaname if it has
        not been closed.
        """
        manager = self.sectionManagers.get(metaname)
        if manager:
            return manager.get_latest_section()

    def appendOnOpen(self, sectionName, onOpen):
        self.sectionManagers.onOpen.append(onOpen)

    def storeToClosedSuper(self, superName, superGIndex, metaInfo, toClose):
        logging.getLogger("nomadcore.caching_backend").warn(
            "Dropping section {name} gIndex {gIndex}, as it cannot be added to closed super section {superName} gIndex: {superGIndex}".format(
                name = metaInfo.name,
                gIndex = toClose.gIndex,
                superName = superName,
                superGIndex = superGIndex))

    def startedParsingSession(self, mainFileUri, parserInfo, parsingStatus = None, parsingErrors = None):
        """should be called when the parsing starts, parserInfo should be a valid json dictionary"""
        if self.propagateStartFinishParsing and self.superBackend:
            self.superBackend.startedParsingSession(mainFileUri, parserInfo, parsingStatus, parsingErrors)

    def finishedParsingSession(self, parsingStatus, parsingErrors, mainFileUri = None, parserInfo = None,
                               parsingStats = None):
        """should be called when the parsing finishes"""
        if self.propagateStartFinishParsing and self.superBackend:
            self.superBackend.finishedParsingSession(parsingStatus, parsingErrors, mainFileUri, parserInfo, parsingStats)

    def addMatchTelemetry(self, match_telemetry, gIndex = -1):
        """ should be called for outputting match telemetry data:
        input data, together with capture info """
        if self.superBackend:
            self.superBackend.addMatchTelemetry(match_telemetry, gIndex)

    def metaInfoEnv(self):
        """the metaInfoEnv this parser was optimized for"""
        return self.__metaInfoEnv

    def openSections(self):
        """returns the sections that are still open
        sections are identified by metaName and their gIndex"""
        res = set()
        for manager in self.sectionManagers:
            res.extend(map(lambda x: (manager.metaInfo.name, x), manager.openSections.keys()))
        return res

    def sectionInfo(self, metaName, gIndex):
        """returns information on a section (for debugging purposes)"""
        manager = self.sectionManagers[metaName]
        section = manager.openSections.get(gIndex, None)
        if section:
            return "open section {} gIndex: {}, references: {}".format(metaName, gIndex, section.references)
        else:
            return "closed section {} gIndex: {}".format(metaName, gIndex)

    def openSection(self, metaName):
        """opens a new section and returns its new unique gIndex"""
        manager = self.sectionManagers[metaName]
        if self.superBackend and manager.forwardOpenClose:
            newGIndex = self.superBackend.openSection(metaName)
            manager.openSectionWithGIndex(self, newGIndex)
        else:
            newGIndex = manager.openSection(self)
        return newGIndex

    def openNonOverlappingSection(self, metaName):
        """opens a new section that is expected not to overlap with itself and returns its new unique gIndex"""
        manager = self.sectionManagers[metaName]
        if len(manager.openSections) != 0:
            raise Exception("Section %s was not expected to overlap" % metaName)
        return self.openSection(metaName)

    def openSectionWithGIndex(self, metaName, gIndex):
        """opens a new section where gIndex is generated externally
        gIndex should be unique (no reopening of a closed section)"""
        manager = self.sectionManagers[metaName]
        if self.superBackend and manager.forwardOpenClose:
            self.superBackend.openSectionWithGIndex(metaName, newGIndex)
        manager.openSectionWithGIndex(newGIndex)

    def setSectionInfo(self, metaName, gIndex, references):
        """sets info values of an open section
        references should be a dictionary with the gIndexes of the root sections this section refers to"""
        manager = self.sectionManagers[metaName]
        manager.setSectionInfo(gIndex, references)
        if self.superBackend and manager.forwardOpenClose:
            self.superBackend.setSectionInfo(metaName, gIndex, references)

    def closeSection(self, metaName, gIndex):
        """closes a section
        after this no other value can be added to the section
        metaName is the name of the meta info, gIndex the index of the section"""
        manager = self.sectionManagers[metaName]
        manager.closeSection(self, gIndex)

    def closeNonOverlappingSection(self, metaName):
        """closes a non overlapping section
        after this no other value can be added to the section
        metaName is the name of the meta info, gIndex the index of the section"""
        manager = self.sectionManagers[metaName]
        openGIndexes = list(manager.openSections.keys())
        if len(openGIndexes) != 1:
            if not openGIndexes:
                raise Exception("Call to closeNonOverlapping(%s) with no open section" % metaName)
            else:
                raise Exception("Section %s was not supposed to overlap, found %s open when closing" % (metaName, openGIndexes))
        manager.closeSection(self, openGIndexes[0])

    def addValue(self, metaName, value, gIndex = -1):
        """adds a json value corresponding to metaName
        the value is added to the section the meta info metaName is in
        a gIndex of -1 means the latest section
        """
        dataManager = self.dataManagers[metaName]
        cachingLevel = dataManager.cachingLevel
        if cachingLevel == CachingLevel.Forward or cachingLevel == CachingLevel.ForwardAndCache:
            if self.superBackend:
                self.superBackend.addValue(metaName, value, gIndex)
        if cachingLevel == CachingLevel.ForwardAndCache or cachingLevel == CachingLevel.Cache or cachingLevel == CachingLevel.PreOpenedCache:
            dataManager.superSectionManager.addValue(dataManager.metaInfo, value, gIndex)

    def addRealValue(self, metaName, value, gIndex = -1, unit=None):
        """Adds a floating point value corresponding to metaName The value is
        added to the section the meta info metaName is in A gIndex of -1 means
        the latest section.
        """
        if unit is not None:
            value = self.convert_unit(metaName, value, unit)

        dataManager = self.dataManagers[metaName]
        cachingLevel = dataManager.cachingLevel
        if cachingLevel == CachingLevel.Forward or cachingLevel == CachingLevel.ForwardAndCache:
            if self.superBackend:
                self.superBackend.addRealValue(metaName, value, gIndex)
        if cachingLevel == CachingLevel.ForwardAndCache or cachingLevel == CachingLevel.Cache or cachingLevel == CachingLevel.PreOpenedCache:
            dataManager.superSectionManager.addValue(dataManager.metaInfo, value, gIndex)

    def addArray(self, metaName, shape, gIndex = -1):
        """adds a new array value of the given size corresponding to metaName
        the value is added to the section the meta info metaName is in
        a gIndex of -1 means the latest section
        the array is unitialized"""
        dataManager = self.dataManagers[metaName]
        cachingLevel = dataManager.cachingLevel
        if cachingLevel == CachingLevel.Forward or cachingLevel == CachingLevel.ForwardAndCache:
            if self.superBackend:
                self.superBackend.addArray(metaName, shape, gIndex)
        if cachingLevel == CachingLevel.ForwardAndCache or cachingLevel == CachingLevel.Cache or cachingLevel == CachingLevel.PreOpenedCache:
            dataManager.superSectionManager.addArray(dataManager.metaInfo, shape, gIndex)

    def setArrayValues(self, metaName, values, offset = None, gIndex = -1, unit=None):
        """Adds values to the last array added, array must be a numpy array
        """
        if unit is not None:
            values = self.convert_unit(metaName, values, unit)

        dataManager = self.dataManagers[metaName]
        cachingLevel = dataManager.cachingLevel
        if cachingLevel == CachingLevel.Forward or cachingLevel == CachingLevel.ForwardAndCache or cachingLevel == CachingLevel.PreOpenedCache:
            if self.superBackend:
                self.superBackend.setArrayValues(metaName, values, offset, gIndex)
        if cachingLevel == CachingLevel.ForwardAndCache or cachingLevel == CachingLevel.Cache or cachingLevel == CachingLevel.PreOpenedCache:
            dataManager.superSectionManager.setArrayValues(dataManager.metaInfo, values, offset, gIndex)

    def addArrayValues(self, metaName, values, gIndex = -1, unit=None):
        """Adds an array value with the given array values.  Values must be a
        numpy array.
        """
        if unit is not None:
            values = self.convert_unit(metaName, values, unit)

        dataManager = self.dataManagers[metaName]
        cachingLevel = dataManager.cachingLevel
        if cachingLevel == CachingLevel.Forward or cachingLevel == CachingLevel.ForwardAndCache:
            if self.superBackend:
                self.superBackend.addArrayValues(metaName, values, gIndex=gIndex)
        if cachingLevel == CachingLevel.ForwardAndCache or cachingLevel == CachingLevel.Cache or cachingLevel == CachingLevel.PreOpenedCache:
            dataManager.superSectionManager.addArrayValues(dataManager.metaInfo, values, gIndex=gIndex)

    def convertScalarStringValue(self, metaName, strValue):
        """converts a scalar string value of the given meta info to a python value"""
        metaInfo = self.metaInfoEnv().infoKindEl(metaName)
        dtypeStr = metaInfo.dtypeStr
        return parser_backend.valueForStrValue(strValue, dtypeStr)

    def arrayForMetaInfo(self, metaName, shape):
        """Returns an array with the correct type for the given meta info, and the given shape"""
        metaInfo = self.metaInfoEnv().infoKindEl(metaName)
        dtypeStr = metaInfo.dtypeStr
        return np.zeros(shape, dtype = parser_backend.numpyDtypeForDtypeStr(dtypeStr))

    def convert_unit(self, metaName, value, unit):
        """Try to perform a unit conversion.
        """
        # If a unit has been specified in the metainfo, check that the
        # conversion can be done
        metaInfo = self.metaInfoEnv().infoKindEl(metaName)
        metainfo_unit = metaInfo.units
        if metainfo_unit:
            metainfo_dim = unit_conversion.ureg(metainfo_unit).dimensionality
            source_dim = unit_conversion.ureg(unit).dimensionality
            if metainfo_dim != source_dim:
                raise Exception("Unit conversion error with metainfo '{}'. Could not convert the given unit '{}' to target unit '{}' specified in the metainfo.".format(metaName, unit, metainfo_unit))

        # If there is a metaname specific unit request, convert to it
        target_unit = None
        if self.metainfo_units is not None:
            target_unit = self.metainfo_units.get(metaName)

        # If there is a dimension specific unit request, convert to it
        if target_unit is None:
            if self.default_units is not None:
                source_dim = unit_conversion.ureg(unit).dimensionality
                map_unit = self.default_units.get(str(source_dim))
                if map_unit:
                    target_unit = map_unit

        converted = unit_conversion.convert_unit(value, unit, target_unit)
        if converted is None:
            raise Exception("ActiveBackend cannot convert to unit '{}'".format(unit))
        return converted

    def pwarn(self, msg):
        """Writes a parser warning message"""
        self.addValue("parsing_message_warning_run", msg)
