#!/usr/bin/env python
# Authors: Franz Knuth & Fawzi Mohamed
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import range
from builtins import object
import json
import numpy as np
import sys, re

def numpyDtypeForDtypeStr(dtypeStr):
    """returns the numpy dtype given the meta info dtype"""
    if dtypeStr == "f" or dtypeStr == "f64":
        return np.float64
    elif dtypeStr == "f32":
        return np.float32
    elif dtypeStr == "i" or dtypeStr == "i32":
        return np.int32
    elif dtypeStr == "i64":
        return np.int64
    elif dtypeStr == "b":
        return bool
    elif dtypeStr == "C":
        return object
    elif dtypeStr == "D":
        return object
    elif dtypeStr == "B":
        return object

def valueForStrValue(strValue, dtypeStr):
    """adds a value translating it from a string"""
    try:
        if dtypeStr[0] == "f":
            return float(strValue.replace("d","e").replace("D", "e"))
        elif dtypeStr[0] == "i":
            return int(strValue)
        elif dtypeStr[0] == "b":
            return (re.match("\s*(?:[nN][oO]?|\.?[fF](?:[aA][lL][sS][eE])?\.?|0)\s*$", strValue) is None)
        elif dtypeStr[0] == "B":
            return strValue # assumed to be base64 encoded
        elif dtypeStr[0] == "C":
            return strValue
        elif dtypeStr[0] == "D":
            return json.loads(strValue)
        else:
            raise Exception("unexpected dtypeStr %s" % (dtypeStr))
    except Exception as e:
        raise Exception("Error when converting %r to dtypeStr %r" % (strValue, dtypeStr), e)

class JsonParseEventsWriterBackend(object):
    """Simple backend that writes out the parse events in json format"""
    # json content is written to fileOut
    def __init__(self, metaInfoEnv, fileOut = sys.stdout, writeMatchTelemetry=False):
        self.__metaInfoEnv = metaInfoEnv
        self.fileOut = fileOut
        self.__gIndex = -1
        self.__openSections = set()
        self.__writeComma = False
        self.__lastIndex = {}
        self.writeMatchTelemetry = writeMatchTelemetry
        self.stats = {}

    def addStat(self, name):
        self.stats[name] = self.stats.get(name, 0) + 1

    def sendStats(self):
        stats = {"parser":{"name":"fhi-aims", "version": "0.3"},
               "data":self.stats}
        url = 'https://nomad-dev.rz-berlin.mpg.de/parsers/addStat'
        #url = 'http://127.0.0.1:8081/parsers/addStat'
        data = json.dumps(stats, sort_keys=True)
        req = urllib.request.Request(url, data)
        response = urllib.request.urlopen(req)
        the_page = response.read()
        sys.stderr.write("stat sending did answer:" + the_page)

    @staticmethod
    def __numpyEncoder(self, o):
        """new default function for json class so that numpy arrays can be encoded"""
        # check if object is a numpy array
        if isinstance(o, np.ndarray):
            # ensure that we have an array with row-major memory order (C like)
            if not o.flags['C_CONTIGUOUS']:
                o = np.ascontiguousarray(o)
            return o.tolist()
            # see default method in python/json/encoder.py
        elif isinstance(o, set):
            return list(sorted(o))
        else:
            raise TypeError(repr(o) + " is not JSON serializable")

    def __jsonOutput(self, dic):
        """method to define format of json output"""
        if self.__writeComma:
            self.fileOut.write(", ")
        else:
            self.__writeComma = True
        json.dump(dic, self.fileOut, indent = 2, separators = (',', ':'), sort_keys=True) # default = self.__numpyEncoder)

    def startedParsingSession(self, mainFileUri, parserInfo, parserStatus = None, parserErrors = None):
        """should be called when the parsing starts, parserInfo should be a valid json dictionary"""
        self.fileOut.write("{\n  \"type\": \"nomad_parse_events_1_0\"")
        self.sessionMainFileUri = mainFileUri
        self.sessionParserInfo = parserInfo
        self.sessionParserStatus = parserStatus
        self.sessionParserErrors = parserErrors
        if mainFileUri is not None:
            self.fileOut.write(",\n  \"mainFileUri\": " + json.dumps(mainFileUri, sort_keys=True))
        if parserInfo is not None:
            self.fileOut.write(",\n  \"parserInfo\": " + json.dumps(parserInfo, indent = 2, separators = (',', ':'), sort_keys=True))
        if parserStatus is not None:
            self.fileOut.write(",\n  \"parserStatus\": " + json.dumps(parserStatus, indent = 2, separators = (',', ':'), sort_keys=True))
        if parserErrors is not None:
            self.fileOut.write(",\n  \"parserStatus\": " + json.dumps(parserErrors, indent = 2, separators = (',', ':'), sort_keys=True))
        self.fileOut.write(""",
  "events": [""")

    def finishedParsingSession(self, parserStatus, parserErrors, mainFileUri = None, parserInfo = None,
                               parsingStats = None):
        """should be called when the parsing finishes"""
        self.fileOut.write("]")
        if mainFileUri is not None and self.sessionMainFileUri is None:
            self.fileOut.write(",\n  \"mainFileUri\": " + json.dumps(mainFileUri, sort_keys=True))
        if parserInfo is not None and self.sessionParserInfo is None:
            self.fileOut.write(",\n  \"parserInfo\": " + json.dumps(parserInfo, indent = 2, separators = (',', ':'), sort_keys=True))
        if parserStatus is not None and self.sessionParserStatus is None:
            self.fileOut.write(",\n  \"parserStatus\": " + json.dumps(parserStatus, indent = 2, separators = (',', ':'), sort_keys=True))
        if parserErrors is not None and self.sessionParserErrors is None:
            self.fileOut.write(",\n  \"parserErrors\": " + json.dumps(parserErrors, indent = 2, separators = (',', ':'), sort_keys=True))
        if parsingStats is not None:
            self.fileOut.write(",\n  \"parsingStats\": " + json.dumps(parsingStats, indent = 4, separators = (',', ':'), sort_keys=True))
        self.fileOut.write("""
}""")
        self.fileOut.flush()

    def openContext(self, contextUri):
        self.__jsonOutput({"event":"openContext", "nomadUri":contextUri})

    def closeContext(self, contextUri):
        self.__jsonOutput({"event":"closeContext", "nomadUri":contextUri})

    def metaInfoEnv(self):
        """the metaInfoEnv this parser was optimized for"""
        return self.__metaInfoEnv

    def openSections(self):
        """returns the sections that are still open
        sections are identified by metaName and their gIndex"""
        return self.__openSections

    def sectionInfo(self, metaName, gIndex):
        """returns information on a section (for debugging purposes)"""
        if (metaName,gIndex) in self.__openSections:
            return "section {} gIndex: {} is open".format(metaName, gIndex)
        else:
            return "section {} gIndex: {} is closed".format(metaName, gIndex)

    def openSection(self, metaName):
        """opens a new section and returns its new unique gIndex"""
        newIndex = self.__lastIndex.get(metaName, -1) + 1
        self.openSectionWithGIndex(metaName, newIndex)
        return newIndex

    def openNonOverlappingSection(self, metaName):
        """opens a new non overlapping section"""
        if any(x[0] == metaName for x in self.__openSections):
            raise Exception("Section %s is not supposed to overlap" % metaName)
        return self.openSection(metaName)

    def openSectionWithGIndex(self, metaName, gIndex):
        """opens a new section where gIndex is generated externally
        gIndex should be unique (no reopening of a closed section)"""
        self.__lastIndex[metaName] = gIndex
        self.__openSections.add((metaName, gIndex))
        self.__jsonOutput({"event":"openSection", "metaName":metaName, "gIndex":gIndex})

    def setSectionInfo(self, metaName, gIndex, references):
        """sets info values of an open section
        references should be a dictionary with the gIndexes of the root sections this section refers to"""
        self.__jsonOutput({"event":"setSectionInfo", "metaName":metaName, "gIndex":gIndex, "references":references})

    def closeSection(self, metaName, gIndex):
        """closes a section
        after this no other value can be added to the section
        metaName is the name of the meta info, gIndex the index of the section"""
        if (metaName, gIndex) in self.__openSections:
            self.__openSections.remove((metaName, gIndex))
            self.__jsonOutput({"event":"closeSection", "metaName":metaName, "gIndex":gIndex})
        # raise exeption if section is not open
        else:
            raise Exception("There is no open section with metaName %s and gIndex %d" % (metaName, gIndex))

    def closeNonOverlappingSection(self, metaName):
        """closes a non overlapping section"""
        openGIndexes = [x for x in self.__openSections if x[0] == metaName]
        if len(openGIndexes) != 1:
            if not openGIndexes:
                raise Exception("Call to closeNonOverlapping(%s) with no open section" % metaName)
            else:
                raise Exception("Section %s was not supposed to overlap, found %s open when closing" % (metaName, openGIndexes))
        self.closeSection(metaName, openGIndexes[0][1])

    def addValue(self, metaName, value, gIndex = -1):
        """adds a json value corresponding to metaName
        the value is added to the section the meta info metaName is in
        a gIndex of -1 means the latest section"""
        self.__jsonOutput({"event":"addValue", "metaName":metaName, "gIndex":gIndex, "value":value})

    def addRealValue(self, metaName, value, gIndex = -1):
        """adds a floating point value corresponding to metaName
        The value is added to the section the meta info metaName is in
        A gIndex of -1 means the latest section"""
        self.__jsonOutput({"event":"addRealValue", "metaName":metaName, "gIndex":gIndex, "value":value})

    def addArray(self, metaName, shape, gIndex = -1):
        """adds a new array value of the given size corresponding to metaName
        the value is added to the section the meta info metaName is in
        a gIndex of -1 means the latest section
        the array is unitialized"""
        self.__jsonOutput({"event":"addArray", "metaName":metaName, "gIndex":gIndex, "shape":shape})

    def setArrayValues(self, metaName, values, offset = None, gIndex = -1):
        """adds values to the last array added, array must be a numpy array"""
        res = {
            "event":"setArrayValues",
            "metaName":metaName,
            "gIndex":gIndex,
            "valuesShape":values.shape,
            "flatValues": values.flatten().tolist()
        }
        if offset:
            res["offset"] = offset
        self.__jsonOutput(res)

    def addArrayValues(self, metaName, values, gIndex = -1):
        """adds an array value with the given array values.
        values must be a numpy array"""
        res = {
            "event":"addArrayValues",
            "metaName":metaName,
            "gIndex":gIndex,
            "valuesShape":values.shape,
            "flatValues": values.flatten().tolist()
        }
        self.__jsonOutput(res)

    def addMatchTelemetry(self, match_telemetry, gIndex = -1):
        if not self.writeMatchTelemetry:
            return
        res = {
            'event': "matchTelemetry",
            'gIndex': gIndex,
            'fInLine': match_telemetry['fInLine'],
            'fInLineNr': match_telemetry['fInLineNr'],
            'matcherName': match_telemetry['matcherName'],
            'matchFlags': match_telemetry['matchFlags'],
            'matchSpansFlat': match_telemetry['matchSpansFlat'],
            'matcherGroup': match_telemetry['matcherGroup'],
        }
        self.__jsonOutput(res)

    def convertScalarStringValue(self, metaName, strValue):
        """converts a scalar string value of the given meta info to a python value"""
        metaInfo = self.metaInfoEnv().infoKindEl(metaName)
        dtypeStr = metaInfo.dtypeStr
        return valueForStrValue(strValue, dtypeStr)

    def arrayForMetaInfo(self, metaName, shape):
        """Returns an array with the correct type for the given meta info, and the given shape"""
        metaInfo = self.metaInfoEnv().infoKindEl(metaName)
        dtypeStr = metaInfo.dtypeStr
        return numpy.zeros(shape, dtype = numpyDtypeForDtypeStr(dtypeStr))

    def pwarn(self, msg):
        """Writes a parser warning message"""
        self.addValue("parsing_message_warning_run", msg)

# testing
if __name__ == '__main__':
    parser = JsonParseEventsWriterBackend(None)

    array = np.array(list(range(3**3)), dtype = 'int64')
    array = np.reshape(array, (3, 3, 3))

    gIndex = parser.openSection('test')
    print("Open sections: " + str(parser.openSections()))
    parser.setSectionInfo('single_run', gIndex, ('main_section', 5))
    parser.addRealValue('energy', 5.0)
    parser.addArrayValues('array', array)
    parser.closeSection('test', gIndex)
    print("Open sections: " + str(parser.openSections()))
