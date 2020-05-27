import setup_paths
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion.unit_conversion import convert_unit
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.metainfo_storage.MetaInfoStorage import COMMON_META_INFO_PATH, PUBLIC_META_INFO_PATH
from nomadcore.metainfo_storage import MetaInfoStorage as mStore
from nomadcore.metainfo_storage.MetaInfoStorage import strcleaner, strisinstance, literal_eval 
from nomadcore.smart_parser.SmartParserDictionary import getDict_MetaStrInDict, getList_MetaStrInDict, get_unitDict
from contextlib import contextmanager
import numpy as np
import json
import os
import re
import io

@contextmanager
def open_section(parser, name):
    gid = parser.openSection(name)
    yield gid
    parser.closeSection(name, gid)

def conv_str(s):
    try:
        return str(s)
    except ValueError:
        return ''

def conv_int(n, default=None):
    try:
        return int(n)
    except ValueError:
        if default:
            return default
        else:
            return None

def conv_float(n, default=None):
    try:
        return float(n)
    except ValueError:
        if default:
            return default
        else:
            return None

def get_metaInfo(self):
    metaInfoEnv, warnings = loadJsonFile(filePath=self.META_INFO_PATH, 
                                         dependencyLoader=None, 
                                         extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS, 
                                         uri=None)
    metaInfoEnv = set_section_metaInfoEnv(self.PARSERTAG, metaInfoEnv, self.sectionDict)
    metaInfoEnv = setDict_metaInfoEnv(metaInfoEnv, self.metaDicts)
    return metaInfoEnv

def set_section_metaInfoEnv(parsertag, infoEnv, sectionDict):
    """Modifies meta info data.

    Args:
        metaInfoEnv: meta info environment json type data.

    Returns:
        metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
    """
    for newName, nameDict in sectionDict.items():
        newName = newName.lower().replace(" ", "").replace("-", "")
        if parsertag+'_%s_%s' % (nameDict["metaNameTag"], newName) not in infoEnv.infoKinds:
            infoEnv.addInfoKindEl(InfoKindEl(
                description='auto generated section meta info data',
                name=parsertag+'_%s_%s' % (nameDict["metaNameTag"], newName),
                kindStr=nameDict["listTypStr"],
                repeats=nameDict["repeatingSection"],
                superNames=nameDict["supraNames"])) 

    return infoEnv

def setDict_metaInfoEnv(infoEnv, metaDicts):
    """Modifies meta info data.

    Args:
        metaInfoEnv: meta info environment json type data.
        nameDict: dictionary for name info and data.

    Returns:
        metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
    """
    for nameList, nameDict in metaDicts.items():
        for keyName in nameDict.keys():
            if '%s' % (keyName) not in infoEnv.infoKinds:
                infoEnv.addInfoKindEl(InfoKindEl(
                    name='%s' % (keyName),
                    description='auto generated meta info data',
                    dtypeStr=nameDict[keyName].metaInfoType,
                    shape=[],
                    superNames=nameDict[keyName].activeSections)) 

    return infoEnv

def set_metaInfoEnv(parsertag, infoEnv, metaNameTag, newList, listTypStr, supraNames):
    """Modifies meta info data.

    Args:
        metaInfoEnv: meta info environment json type data.

    Returns:
        metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
    """
    for newName in newList:
        newName = newName.lower().replace(" ", "").replace("-", "")
        if 'x_'+parsertag+'_%s_%s' % (metaNameTag, newName) not in infoEnv.infoKinds:
            infoEnv.addInfoKindEl(InfoKindEl(
                name='x_'+parsertag+'_%s_%s' % (metaNameTag, newName),
                description='auto generated meta info data',
                dtypeStr=listTypStr,
                shape=[],
                superNames=supraNames)) 

    return infoEnv

def namelist_matcher(cText, key, matchWith):
    matchThisParsy = None
    if 'EOL' in matchWith:
        matchThisParsy = re.compile(r"(?:%s)\s*(?:\s|=|:)\s*(?:'|\")?"
                                     "(?P<%s>.*)(?:'|\")?\s*,?"
                                     % (cText, key))
    elif 'UD' in matchWith:
        delimeter = matchWith.replace('UD', '')
        matchThisParsy = re.compile(r"(?:%s)\s*(?:\s|=|:)\s*(?:'|\")?"
                                     "(?P<%s>[\-+0-9.a-zA-Z:]+)"
                                     "(?:'|\")?\s*[%s]" 
                                     % (cText, key, delimeter))
    elif 'BOL' in matchWith:
        matchThisParsy = re.compile(r"(?:'|\")?(?P<%s>.*)(?:'|\")?"
                                     "\s*(?:%s)\s*"
                                     % (key, cText))
    elif 'PW' in matchWith:
        matchThisParsy = re.compile(r"(?:'|\")?(?P<%s>[\-+0-9.a-zA-Z:]+)"
                                     "(?:'|\")?\s*(?:%s)\s*"
                                     % (key, cText))
    elif 'MW' in matchWith:
        matchThisParsy = re.compile(r"(?P<%s>[\S]*(?=%s)[\S]*)" 
                                     % (key, cText))
    elif 'FD' in matchWith:
        delimeter = matchWith.replace('FD', '')
        matchThisParsy = re.compile(r"[%s]\s*(?:'|\")?"
                                     "(?P<%s>[\-+0-9.a-zA-Z:]+)"
                                     "(?:'|\")?\s*(?:%s)\s*"
                                     % (delimeter, key, cText))
    else: # Default matchWith 'NW'
        matchThisParsy = re.compile(r"(?:%s)\s*(?:\s|=|:)\s*(?:'|\"|{|\[|\()?"
                                     "(?P<%s>[\S]+)"
                                     "(?:'|\"|}|\]|\))?\s*" 
                                     % (cText, key))
    return matchThisParsy

class ParserBase(object):
    """Base class for SmartParsers"""
    def __init__(self,cachingLevelForMetaName=None, coverageIgnoreList=None,
                 re_program_name=None, parsertag=None, metainfopath=None,
                 parserinfodef=None, recorderOn=False):
        self.PARSERTAG = parsertag
        self.metaStorage = mStore.Container('section_run')
        self.metaStorageRestrict = mStore.Container('section_restricted_uri')
        exclude_dict = { 
            'section_run' : [
            'section_processor_info', 
            'section_processor_log', 
            'section_springer_material',
            'section_repository_info'
            ]}
        self.META_INFO_PATH = metainfopath
        PARSER_INFO_DEFAULT = parserinfodef
        jsonmetadata = mStore.JsonMetaInfo(
                COMMON_META_INFO_PATH, 
                PUBLIC_META_INFO_PATH,
                self.META_INFO_PATH
                )
        self.metaStorage.build(jsonmetadata, 'section_run', exclude_dict)
        self.metaStorageRestrict.build(jsonmetadata, 'section_restricted_uri', exclude_dict)
        self.re_program_name = re_program_name
        #set_Dictionaries(self)
        self.parserInfo = PARSER_INFO_DEFAULT.copy()
        #self.metaInfoEnv = get_metaInfo(self)
        self.coverageIgnoreList = [
            # ignore empty lines
            r"\s*",
            # table separators
            #r"^\s*[=%-]+\s*$",
            #r"^\s*%\s*%\s*$",
        ]
        self.coverageIgnore = None
        self.strcleaner = strcleaner
        self.strisinstance = strisinstance
        self.literal_eval = literal_eval
        self.recordList = None
        if recorderOn:
            self.recordList = io.StringIO()

    def parse(self):
        self.coverageIgnore = re.compile(r"^(?:" + r"|".join(self.coverageIgnoreList) + r")$")
        mainFunction(mainFileDescription=self.mainFileDescription(), 
                     metaInfoEnv=self.metaInfoEnv, 
                     parserInfo=self.parserInfo,
                     cachingLevelForMetaName=self.cachingLevelForMetaName,
                     superContext=self)
        if self.recordList:
            self.recordList.close()

    def peekline(self, parser):
        pos = parser.fIn.fIn.tell()
        line = parser.fIn.fIn.readline()
        parser.fIn.fIn.seek(pos)
        return line

    def peekrecord(self):
        pos = self.recordList.tell()
        line = self.recordList.readline()
        self.recordList.seek(pos)
        return line

    def topology_system_name(self, backend, gIndex):
        """ Function to generate data for system_name
        """
        system_name=None
        supB = backend.superBackend
        if self.topology:
            topoDict = self.topology.topoDict
            if "system_name" in topoDict:
                system_name = topoDict["system_name"]
        if system_name is not None:
            supB.addValue('system_name', str(system_name))    

    def topology_atom_type_and_interactions(self, backend, gIndex):
        """ Function to generate data for atom_to_molecule 
        """
        sO = open_section
        supB = backend.superBackend
        topoDict = None
        if self.topology:
            topoDict = self.topology.topoDict
        if topoDict:
            if "types_list" in topoDict:
                types_list = topoDict["types_list"]
            else:
                types_list = []
            if "atom_type_list" in topoDict:
                atom_type_list = topoDict["atom_type_list"]
            else:
                atom_type_list = []
            if "atom_label_list" in topoDict:
                atomlabelList = topoDict["atom_label_list"]
            else:
                atomlabelList = []
            if "name_list" in topoDict:
                name_list = topoDict["name_list"]
            else:
                name_list = []
            if "atom_mass_list" in topoDict:
                massesList = topoDict["atom_mass_list"]
            else:
                massesList = []
            if "element_list" in topoDict:
                elementList = topoDict["element_list"]
            else:
                elementList = []
            if "atom_element_list" in topoDict:
                atomelementList = topoDict["atom_element_list"]
            else:
                atomelementList = []
            if "atom_radius_list" in topoDict:
                radiusList = topoDict["atom_radius_list"]
            else:
                radiusList = []
            if "atom_charge_list" in topoDict:
                chargeList = topoDict["atom_charge_list"]
            else:
                chargeList = []
            if "interactions_dict" in topoDict:
                interDict = topoDict["interactions_dict"]
            else:
                interDict = []
            if "interactions_type_dict" in topoDict:
                interTypeDict = topoDict["interactions_type_dict"]
            else:
                interTypeDict = []
            for elm in range(len(types_list)):
                with sO(supB, 'section_atom_type'):
                    # Atom Type NOMAD Name 
                    # Q: Why different than name or type? 
                    # A: Because atom name in topology is used to define a specific atom 
                    #    with properties may be same or different with the other name-sharers. 
                    #    On the other hand, atom type defines a type that mostly fits 
                    #    with the expectations of a Force Field parameter set where a 
                    #    specific atom-atom pair or atom-mutliatom set have an interaction. 
                    #    A unique identification of an atom may be in both worlds not solely 
                    #    atom type or atom name. While there are cases where atom name, 
                    #    atom element symbol and even atom type are same. Such cases happen 
                    #    quiet often in ab-inito calculations.
                    #    But this is not the case that we are talking about here in MD 
                    #    calculations since what we expect from a NOMAD meta info that 
                    #    it should have a lossless translation of a topology not more.
                    #
                    supB.addValue('atom_type_name', str(types_list[elm]) + ':' + str(elm+1))    
                    # Atom Name (Topology Name)
                    supB.addValue(self.PARSERTAG + '_atom_name', name_list[elm]) 
                    # Atom Type Name (Force Field Name)
                    supB.addValue(self.PARSERTAG + '_atom_type', atom_type_list[elm]) 
                    # Atomic element/symbol
                    supB.addValue(self.PARSERTAG + '_atom_element', elementList[elm]) 
                    # Atomic mass
                    if massesList[elm] is not None:
                        supB.addValue('atom_type_mass', massesList[elm])             
                    if radiusList[elm]:
                        # Atomic van der Waals radius
                        supB.addValue(self.PARSERTAG + '_atom_radius', radiusList[elm])   
                    if chargeList[elm]:
                        # Atomic charge
                        supB.addValue('atom_type_charge', chargeList[elm])
#                    if atomlabelList[elm]:
#                        supB.addArrayValues(self.PARSERTAG + '_atom_to_atom_type_ref', 
#                                np.asarray(atomlabelList[elm]))   
                    pass

            with sO(supB, self.PARSERTAG + '_section_atom_to_atom_type_ref'):
                if atomlabelList:
                    supB.addArrayValues(self.PARSERTAG + '_atom_to_atom_type_ref', 
                            np.asarray(atomlabelList))   

            for inum in range(len(interDict.keys())):
                with sO(supB, 'section_interaction'):
                    # atom indexes of bound pairs for a specific atom type
                    supB.addArrayValues('interaction_atoms', np.asarray(interDict[inum]))     
                    # number of bonds for this type
                    supB.addValue('number_of_interactions', len(interDict[inum]))             
                    # number of atoms involved (2 for bonds, 3 for angles ,and so on)
                    supB.addValue('number_of_atoms_per_interaction', len(interDict[inum][0])) 
                
                    #if bondFunctional:
                    #    self.superP.addValue('interaction_kind', bondFunctional)  # functional form of the interaction
                    # this points to the relative section_atom_type
                    supB.addArrayValues(self.PARSERTAG + '_interaction_atom_to_atom_type_ref', 
                            np.asarray(interTypeDict[inum]))  
                    # interaction parameters for the functional
                    #self.superP.addValue('interaction_parameters', bondParameters)  


#            for i in range(len(moleculeTypeInfo)):
#
#                with o(p, 'section_molecule_type'):
#                    # gindex = 0
#                    p.addValue('molecule_type_name', 'molecule'+'_'+str(moleculeTypeInfo[i][0]))
#                    p.addValue('number_of_atoms_in_molecule', len(moleculeTypeInfo[i][1]))
#
#                    p.addArrayValues('atom_in_molecule_to_atom_type_ref', np.asarray([x-1 for x in moleculeTypeInfo[i][1]]))
#

#
#    mol_to_mol_bond = []
#    mol_to_mol_bond_num = []
#    atom_in_mol_bond = []
#    atom_in_mol_bond_num = []
#    index=0
#    print(len([res for res in mytopology.bonds]))
#    for res in mytopology.bonds:
#        molname1, molatom1 = str(res[0]).split('-')
#        molname2, molatom2 = str(res[1]).split('-')
#        if molname1 in molname2: 
#            atom_in_mol_bond.append(res)
#            atom_in_mol_bond_num.append(bonds[index])
#        else:
#            mol_to_mol_bond.append(res)
#            mol_to_mol_bond_num.append(bonds[index])
#        index += 1
#    print(mol_to_mol_bond)
#    print(np.array(mol_to_mol_bond_num))

    def adHoc_program_name(self, parser):
        if self.re_program_name is not None:
            if not self.re_program_name.match(
                    parser.lastMatch[self.PARSERTAG+'_program_name']):
                raise Exception(
                    "mainFile program name was: %s, unsuited for %s" % (
                        parser.lastMatch[self.PARSERTAG+'_program_name'],
                        type(self).__name__))

    def dictionary_parser(self, parser, stopOnMatchStr, quitOnMatchStr, metaNameStart, matchNameList, 
            matchNameDict, updateMatchDict, onlyCaseSensitive, stopOnFirstLine, parserOptions, 
            entryline=None, parserID=None, parsername=None, globalDict=None, localDict=None, parent=None):
        record = None
        replay = None
        replayCount = None
        if parent is None:
            parent = 0
        rank = parent + 1
        if globalDict is not None:
            record = globalDict["record"]
            replay = globalDict["replay"]
            replayCount = globalDict["replayCount"]
        if record is None:
            record = False
        if replayCount is None:
            replayCount = 0
        if entryline is None:
            lastLine = parser.fIn.fInLine
        else:
            lastLine = entryline
        parserDict = {
                "firstLine"   : 0,
                "storedLines" : '',
                "numStoredLines" : 0,
                "parserID" : parserID,
                "parent" : parent,
                "rank" : rank
                }
        parserDict.update({"firstLine" : parserDict["firstLine"] + 1})
        parserSuccess = False
        dictionary = None
        dicttype = True
        readdict = True
        keyMapper = None
        updatefunc = None
        updateattrs = None
        updateconvert = None
        cntrlattrs = None
        preprocess = None
        postprocess = None
        parsercntrlattr = None
        parsercntrlin = None
        lookupdict = None
        lookupvals = None
        rnewline = False
        checkdict = {}
        if "readline" in parserOptions.keys():
            if parserOptions["readline"] == True:
                rnewline = True
        if "dictionary" in parserOptions.keys():
            if isinstance(parserOptions["dictionary"], dict):
                dictionary = parserOptions["dictionary"] 
            if isinstance(parserOptions["dictionary"], str):
                dictionary = getattr(self,parserOptions["dictionary"]) 
            else:
                dictionary = None
        if "dicttype" in parserOptions.keys():
            dicttype = True if "smartparser" in parserOptions["dicttype"] else False
        if "readwritedict" in parserOptions.keys():
            if parserOptions["readwritedict"] == "write":
                readdict = False
        if "keyMapper" in parserOptions.keys():
            keyMapper = parserOptions["keyMapper"] if isinstance(parserOptions["keyMapper"], dict) else None
        if "updatefunc" in parserOptions.keys():
            if isinstance(parserOptions["updatefunc"], str):
                if 'max' in parserOptions["updatefunc"]:
                    updatefunc = np.amax
                    updateconvert = float
                if 'min' in parserOptions["updatefunc"]:
                    updatefunc = np.amin
                    updateconvert = float
            else:
                updatefunc = parserOptions["updatefunc"]
        if "updateattrs" in parserOptions:
            updateattrs = parserOptions["updateattrs"]
        if "controlattrs" in parserOptions:
            cntrlattrs = parserOptions["controlattrs"]
            if cntrlattrs:
                for cntrla in cntrlattrs:
                    if hasattr(self,cntrla):
                        checkdict.update({cntrla:getattr(self,cntrla)})
            #print("PRINTING 0 checkdict:",checkdict)
        if "preprocess" in parserOptions.keys():
            preprocess = parserOptions["preprocess"]
        if "postprocess" in parserOptions.keys():
            postprocess = parserOptions["postprocess"]
        if "lookupdict" in parserOptions:
            lookupdict = getattr(self,parserOptions["lookupdict"])
        if "parsercntrlattr" in parserOptions:
            parsercntrlattr = getattr(self,parserOptions["parsercntrlattr"])
        if "parsercntrlin" in parserOptions:
            if lookupdict:
                parsercntrlin = parserOptions["parsercntrlin"]
                if parsercntrlin in lookupdict:
                    lookupvals = lookupdict[parsercntrlin]
            else:
                if isinstance(parserOptions["parsercntrlin"], str):
                    parsercntrlin = getattr(self,parserOptions["parsercntrlin"])
                else:
                    parsercntrlin = parserOptions["parsercntrlin"]

        mNameDict = getattr(self, matchNameDict)

        if dictionary is not None:
            if dicttype:
                dictionaryStr = getDict_MetaStrInDict(dictionary)
            else:
                dictionaryStr = dictionary
            anythingmatched = False
            # Search for dictionary keys
            for cName, key in getDict_MetaStrInDict(mNameDict).items():
                cNewNameList = []
                # Convert with filter
                if keyMapper is not None:
                    for k, v in keyMapper.items():
                        cNameBack = cName
                        cNewNameList.append(cNameBack.replace(k, v))
                else:
                    cNewNameList.append(cName)
                for cNewName in cNewNameList:
                    if cNewName in dictionaryStr:
                        # if parser control is not given
                        # default is first element of the input
                        elementid = 0
                        if lookupvals is not None:
                            if parsercntrlattr in lookupvals:
                                elementid = lookupvals.index(parsercntrlattr)
                        if readdict:
                            if dicttype:
                                val = dictionary[dictionaryStr[cNewName]].value
                            else:
                                val = dictionary[cNewName]
                            if isinstance(val, (list,tuple,np.ndarray)):
                                v = val[elementid]
                            else:
                                v = val
                            if preprocess is not None:
                                v = preprocess(cName,v)
                            if updatefunc is not None:
                                if updateconvert is None:
                                    v=updatefunc(v)
                                else:
                                    vallist=[]
                                    if mNameDict[key].value is not None:
                                        vallist.append(updateconvert(mNameDict[key].value))
                                    if v is not None:
                                        vallist.append(updateconvert(v))
                                    if not vallist:
                                        pass
                                    else:
                                        v=updatefunc(np.asarray(vallist))
                            if postprocess is not None:
                                v = postprocess(cName,v)
                            if key in list(parser.lastMatch.keys()):
                                parser.lastMatch[key]=v
                            else:
                                mNameDict[key].value=v
                                mNameDict[key].activeInfo=True
                            if isinstance(updateattrs,(tuple,list)):
                                for uattr in updateattrs:
                                    if(uattr == cNewName and hasattr(self,cNewName)):
                                        setattr(self,uattr,v)
                                        #print("PRINTING: self update2:",uattr,v)
                        else:
                            val = mNameDict[key].value
                            #print("PRINTING: dict key val:",key,val)
                            if isinstance(val, (list,tuple,np.ndarray)):
                                v = val[elementid]
                            else:
                                v = val
                            if preprocess is not None:
                                #v = preprocess(v)
                                v = preprocess(cName,v)
                            if updatefunc is not None:
                                vallist=[]
                                if dicttype:
                                    if dictionary[dictionaryStr[cNewName]].value is not None:
                                        vallist.append(updateconvert(
                                            dictionary[dictionaryStr[cNewName]].value
                                            ))
                                else:
                                    if dictionary[cNewName] is not None:
                                        vallist.append(updateconvert(
                                            dictionary[cNewName]
                                            ))
                                if v is not None:
                                    vallist.append(updateconvert(v))
                                if not vallist:
                                    pass
                                else:
                                    #print("PRINTING: vallist:",vallist)
                                    v=updatefunc(np.asarray(vallist))
                            if postprocess is not None:
                                v = postprocess(cName,v)
                            if dicttype:
                                dictionary[dictionaryStr[cNewName]].value=v
                                #print("PRINTING: update key val:",dictionaryStr[cNewName],v)
                            else:
                                dictionary[cNewName]=v
                                #print("PRINTING: update key val:",cNewName,v)
                        anythingmatched = True
                        if readdict:
                            pass
                        else:
                            setattr(self,parserOptions["dictionary"],dictionary)
                            if isinstance(updateattrs,(tuple,list)):
                                for dk,dv in dictionary.items():
                                    for uattr in updateattrs:
                                        if(uattr == dk and hasattr(self,uattr)):
                                            setattr(self,uattr,dv)
                                            #print("PRINTING: self uupdate2:",uattr,dv)

            cntrlcheck = False
            if cntrlattrs:
                for cntrla in cntrlattrs:
                    if hasattr(self,cntrla):
                        #print("PRINTING getattr checkdict:",cntrla,getattr(self,cntrla),checkdict[cntrla])
                        if getattr(self,cntrla)==checkdict[cntrla]:
                            cntrlcheck = True
                            break

            if cntrlcheck:
                # We have matched keywords to update at sections
                # if the active sections are defined
                cntrlDict = None
                cntrlsec = None
                cntrlsave = None
                if "controldict" in parserOptions:
                    cntrlDict = getattr(self,parserOptions["controldict"])
                if "controlsections" in parserOptions:
                    cntrlsec = parserOptions["controlsections"]
                if "controlsave" in parserOptions:
                    if cntrlDict:
                        if parserOptions["controlsave"] in cntrlDict:
                            secDict = {}
                            cntrlsave = cntrlDict[parserOptions["controlsave"]]
                            if isinstance(cntrlsave,dict):
                                secDict = cntrlsave
                            if cntrlsec is not None:
                                for sec in cntrlsec:
                                    secDict.update({sec:True})
                            else:
                                addname = "table_store"
                                if parserDict["parserID"] is not None:
                                    addname = addname + str(parserDict["parserID"])
                                secDict.update({addname:True})
                            cntrlDict.update({parserOptions["controlsave"]:secDict})
                    else:
                        secDict = {}
                        cntrlsave = getattr(self,parserOptions["controlsave"])
                        if isinstance(cntrlsave,dict):
                            secDict = cntrlsave
                        if cntrlsec is not None:
                            for sec in cntrlsec:
                                secDict.update({sec:True})
                        else:
                            addname = "table_store"
                            if parserDict["parserID"] is not None:
                                addname = addname + str(parserDict["parserID"])
                            secDict.update({addname:True})
                        setattr(self,parserOptions["controlsave"],secDict)
                    setattr(self,parserOptions["controldict"],cntrlDict)
                parserSuccess = True

        if(mNameDict is not None and updateMatchDict):
            setattr(self, matchNameDict, mNameDict)
        #print("PRINTING: dict Md step:",self.MDcurrentstep)
        if rnewline is True:
            # Playing from record?
            if record and replayCount>0:
                lastLine = self.recordList.readline()
            else:
                lastLine = parser.fIn.readline()
                # If not replaying, are we recording?
                if record:
                    self.recordList.write(lastLine)
        return lastLine, parserSuccess, globalDict, localDict

    def readline_control_parser(self, parser, stopOnMatchStr, quitOnMatchStr, metaNameStart, matchNameList, 
            matchNameDict, updateMatchDict, onlyCaseSensitive, stopOnFirstLine, parserOptions, entryline=None, 
            parserID=None, parsername=None, globalDict=None, localDict=None, parent=None):
        record = None
        replay = None
        replayCount = None
        if parent is None:
            parent = 0
        rank = parent + 1
        if globalDict is not None:
            record = globalDict["record"]
            replay = globalDict["replay"]
            replayCount = globalDict["replayCount"]
        if record is None:
            record = False
        if replayCount is None:
            replayCount = 0
        if entryline is None:
            lastLine = parser.fIn.fInLine
        else:
            lastLine = entryline
        parserDict = {
                "firstLine"   : 0,
                "storedLines" : '',
                "numStoredLines" : 0,
                "parserID" : parserID,
                "parent" : parent,
                "rank" : rank
                }
        parserSuccess = False
        peeklineFirst = False
        waitatlineRe = None
        controlattr = None
        controlnextattr = None
        controllast = None
        controlskip = []
        controlin = []
        controlwait = None
        controlcounter = 0
        controldict = None
        lookupdict = None
        stopOnMatchRe = None
        quitOnMatchRe = None
        if stopOnMatchStr is not None:
            stopOnMatchRe = re.compile(stopOnMatchStr)
        if quitOnMatchStr is not None:
            quitOnMatchRe = re.compile(quitOnMatchStr)
        if "waitatlineStr" in parserOptions:
            if parserOptions["waitatlineStr"] is not None:
                waitatlineRe = re.compile(parserOptions["waitatlineStr"])
        if "controldict" in parserOptions:
            controldict = getattr(self,parserOptions["controldict"])
        if "lookupdict" in parserOptions:
            lookupdict = getattr(self,parserOptions["lookupdict"])
        if "controlattr" in parserOptions:
            controlattr = getattr(self,parserOptions["controlattr"])
        if "controlnextattr" in parserOptions:
            controlnextattr = getattr(self,parserOptions["controlnextattr"])
        if "controllast" in parserOptions:
            controllast = parserOptions["controllast"]
        if "peeklineFirst" in parserOptions:
            peeklineFirst = parserOptions["peeklineFirst"]
        if "controlin" in parserOptions:
            if lookupdict:
                if parserOptions["controlin"] in lookupdict:
                    controlin = lookupdict[parserOptions["controlin"]]
            else:
                if isinstance(parserOptions["controlin"], str):
                    controlin = getattr(self,parserOptions["controlin"])
                else:
                    controlin = parserOptions["controlin"]
            if("controlskip" in parserOptions and 
                controlin is not None):
                controlskip = [
                        controlin[cs] for cs in parserOptions[
                            "controlskip"] if cs < len(controlin)
                        ]
        if "controlwait" in parserOptions:
            if lookupdict:
                if parserOptions["controlwait"] in lookupdict:
                    controlwait = lookupdict[parserOptions["controlwait"]]
            else:
                if isinstance(parserOptions["controlwait"], str):
                    controlwait = getattr(self,parserOptions["controlwait"])
                else:
                    controlwait = parserOptions["controlwait"]
        if "controlcounter" in parserOptions:
            if controldict:
                if parserOptions["controlcounter"] in controldict:
                    controlcounter = controldict[parserOptions["controlcounter"]]
                else:
                    controldict.update({parserOptions["controlcounter"] : controlnextattr})

        continuenewline = False
        if controlattr in controlskip:
            continuenewline = True
            if "controlcounter" in parserOptions:
                if controldict:
                    controldict.update({parserOptions["controlcounter"] : controlnextattr})
        elif controlattr > controlin[-1]:
            continuenewline = True
            if "controlcounter" in parserOptions:
                if controldict:
                    controldict.update({parserOptions["controlcounter"] : controlnextattr})
        else:
            if controlattr == controlcounter:
                continuenewline = True
                if "controlcounter" in parserOptions:
                    if controldict:
                        controldict.update({parserOptions["controlcounter"] : controlnextattr})
            else:
                continuenewline = False
                if waitatlineRe:
                    if waitatlineRe.findall(lastLine):
                        continuenewline = False
                    else:
                        if controlwait is not None:
                            if controlattr in controlwait:
                                continuenewline = True
                        else:
                            continuenewline = True
                else:
                    if controlwait is not None:
                        if controlattr in controlwait:
                            continuenewline = True
                    #else:
                    #    continuenewline = True
                if continuenewline is False:
                    if controlwait is not None:
                        if controlattr in controlwait:
                            continuenewline = True
                    else:
                        continuenewline = True

        #print("PRINTING nextlogsteps, MDnextstep, targetstep:",controlwait,controlnextattr,controlcounter)
        #print("PRINTING continuenewline:",continuenewline)
        
        if continuenewline:
            #print("PRINTING: readline_control lastLine:",lastLine)
            parserSuccess=True
            peekedline=None
            movenewline=True
            if peeklineFirst is True:
                if "peekedline" in localDict:
                    peekedline=localDict["peekedline"]
                # Playing from record?
                if record and replayCount>0:
                    lastLine = self.peekrecord()
                else:
                    lastLine = self.peekline(parser)
                if peekedline is None:
                    movenewline=False
                    localDict.update({"peekedline" : lastLine})
                    #print("PRINTING: 1 peeklineFirst: ",globalDict["peekedline"],movenewline)
                else:
                    if peekedline == lastLine:
                        movenewline=True
                #print("PRINTING: 2 peeklineFirst: ",peekedline,movenewline)
            if movenewline is True:
                if stopOnMatchStr is not None:
                    if stopOnMatchRe.findall(lastLine):
                        #print("PRINTING: movenewline matchRe:",rank,stopOnMatchRe.findall(lastLine),lastLine)
                        movenewline = False
                if quitOnMatchRe is not None:
                    if quitOnMatchRe.findall(lastLine):
                        #print("PRINTING: check_subparser Rank QuitMatched:",quitOnMatchRe.findall(lastLine),lastLine)
                        movenewline = False
            if movenewline is True:
                # Playing from record?
                if record and replayCount>0:
                    lastLine = self.recordList.readline()
                else:
                    lastLine = parser.fIn.readline()
                    # If not replaying, are we recording?
                    if record:
                        self.recordList.write(lastLine)
                localDict.update({"peekedline" : lastLine})
        else:
            # Playing from record?
            if record and replayCount>0:
                lastLine = self.peekrecord()
            else:
                lastLine = self.peekline(parser)
                # No need to record peeklines
            localDict.update({"peekedline" : lastLine})
        parserDict.update({"firstLine" : parserDict["firstLine"] + 1})

        if controldict:
            setattr(self,parserOptions["controldict"],controldict)
        return lastLine, parserSuccess, globalDict, localDict

    def section_control_parser(self, parser, stopOnMatchStr, quitOnMatchStr, metaNameStart, matchNameList, 
            matchNameDict, updateMatchDict, onlyCaseSensitive, stopOnFirstLine, parserOptions, entryline=None, 
            parserID=None, parsername=None, globalDict=None, localDict=None, parent=None):
        if parent is None:
            parent = 0
        rank = parent + 1
        if globalDict is not None:
            record = globalDict["record"]
            replayCount = globalDict["replayCount"]
        if entryline is None:
            lastLine = parser.fIn.fInLine
        else:
            lastLine = entryline
        parserDict = {
                "firstLine"   : 0,
                "storedLines" : '',
                "numStoredLines" : 0,
                "parserID" : parserID,
                "parent" : parent,
                "rank" : rank
                }
        if parsername is None:
            parsername = "section_control_parser_"+str(parserID)
        sectionopenname = parsername + "_open"
        sectionclosename = parsername + "_close"
        parserSuccess = False
        backend = parser.backend
        sectionname = None
        sectionopen = False
        sectionopenattr = None
        sectionopenin = None
        sectionopenrecord = None
        sectionclose = False
        sectioncloseattr = None
        sectionclosein = None
        sectioncloserecord = None
        lookupdict = None
        activatesection = None
        if "sectionname" in parserOptions:
            sectionname = parserOptions["sectionname"]
        if "lookupdict" in parserOptions:
            lookupdict = getattr(self,parserOptions["lookupdict"])
        if "sectionopen" in parserOptions:
            sectionopen = parserOptions["sectionopen"]
        if "sectionopenattr" in parserOptions:
            sectionopenattr = getattr(self,parserOptions["sectionopenattr"])
        if "sectionopenin" in parserOptions:
            if lookupdict:
                if parserOptions["sectionopenin"] in lookupdict:
                    sectionopenin = lookupdict[parserOptions["sectionopenin"]]
            else:
                if isinstance(parserOptions["sectionopenin"], str):
                    sectionopenin = getattr(self,parserOptions["sectionopenin"])
                else:
                    sectionopenin = parserOptions["sectionopenin"]
            if sectionopenin is not None:
                if lookupdict:
                    if sectionopenname in lookupdict:
                        sectionopenrecord = lookupdict[sectionopenname]
                else:
                    sectionopenrecord = getattr(self,sectionopenname)
                if sectionopenrecord is None:
                    #print("PRINTING: section open record setup")
                    sectionopenrecord = dict()
                    for item in sectionopenin:
                        sectionopenrecord.update({str(item):False})
        if "sectionclose" in parserOptions:
            sectionclose = parserOptions["sectionclose"]
        if "sectioncloseattr" in parserOptions:
            sectioncloseattr = getattr(self,parserOptions["sectioncloseattr"])
        if "sectionclosein" in parserOptions:
            if lookupdict:
                if parserOptions["sectionclosein"] in lookupdict:
                    sectionclosein = lookupdict[parserOptions["sectionclosein"]]
            else:
                if isinstance(parserOptions["sectionclosein"], str):
                    sectionclosein = getattr(self,parserOptions["sectionclosein"])
                else:
                    sectionclosein = parserOptions["sectionclosein"]
            if sectionclosein is not None:
                if lookupdict:
                    if sectionclosename in lookupdict:
                        sectioncloserecord = lookupdict[sectionclosename]
                else:
                    sectioncloserecord = getattr(self,sectionclosename)
                if sectioncloserecord is None:
                    #print("PRINTING: section close record setup")
                    sectioncloserecord = dict()
                    for item in sectionclosein:
                        sectioncloserecord.update({str(item):False})
        if "activatesection" in parserOptions:
            if lookupdict:
                if parserOptions["activatesection"] in lookupdict:
                    activatesection = lookupdict[parserOptions["activatesection"]]

        if sectionname:
            activate = False
            if activatesection is not None:
                if sectionname in activatesection:
                    activate = activatesection[sectionname]
            if activate:
                if(sectionopen is not None and 
                    sectionopenattr is not None and 
                    sectionopenrecord is not None):
                    if str(sectionopenattr) in sectionopenrecord:
                        if sectionopenrecord[str(sectionopenattr)] is False:
                            parserSuccess = True
                            gIndex = backend.openSection(sectionname)
                            self.secGIndexDict.update({sectionname : [1, gIndex]})
                            sectionopenrecord[str(sectionopenattr)] = True
                if(sectionclose is not None and 
                    sectioncloseattr is not None and 
                    sectioncloserecord is not None):
                    #print("PRINTING: section close:",str(sectioncloseattr))
                    if str(sectioncloseattr) in sectioncloserecord:
                        if sectioncloserecord[str(sectioncloseattr)] is False:
                            if self.secGIndexDict[sectionname][0]:
                                parserSuccess = True
                                backend.closeSection(sectionname, self.secGIndexDict[sectionname][1])
                                self.secGIndexDict.update({sectionname : [0, gIndex]})
                                sectioncloserecord[str(sectioncloseattr)] = True

        if lookupdict:
            lookupdict.update({sectionopenname:sectionopenrecord})
            lookupdict.update({sectionclosename:sectioncloserecord})
        else:
            setattr(self,sectionopenname,sectionopenrecord)
            setattr(self,sectionclosename,sectioncloserecord)

        if parserSuccess:
            if activatesection is not None:
                if sectionname in activatesection:
                    activatesection[sectionname] = False
                    lookupdict.update({parserOptions["activatesection"] : activatesection})
                setattr(self,parserOptions["lookupdict"],lookupdict)

        return lastLine, parserSuccess, globalDict, localDict

    def table_store(self, parser, lastLine, stopOnMatchRe, quitOnMatchRe, 
            metaNameStart, matchNameList, matchNameDict, updateMatchDict, 
            onlyCaseSensitive, stopOnFirstLine, parserDict, parserOptions, 
            globalDict=None, localDict=None):
        """
            header = False     if True, the header titles will be read
            wrap = False       if True, table titles and values are assumed to be wrapped
            tablelines = 0     Number of lines that are cycled in the table (Default: 0= process 
                               every line) Ex.: 1= process after 1 line, useful if table is wraped
            lineFilter = None  A pass through (converter) filter before parsing
            headerList = None  if header is not True, values will be assigned to 
                               titles in the headerList in order as seen in line
        """
        record = None
        replay = None
        replayCount = None
        if globalDict is not None:
            record = globalDict["record"]
            replay = globalDict["replay"]
            replayCount = globalDict["replayCount"]
        if record is None:
            record = False
        if replayCount is None:
            replayCount = 0
        if onlyCaseSensitive is None:
            onlyCaseSensitive = True
        #print("PRINTING: table_store lastLine:",lastLine)
        tableStartRe = None
        if "tablestartsat" in parserOptions:
            tableStartRe = re.compile(parserOptions["tablestartsat"])
        tableEndRe = None
        if "tableendsat" in parserOptions:
            tableEndRe = re.compile(parserOptions["tableendsat"])
        stopOnMatch = False
        if isinstance(stopOnMatchRe, str):
            pass
        else:
            if stopOnMatchRe.findall(lastLine):
                stopOnMatch = True
                if parserDict["firstLine"] == 0:
                    if stopOnFirstLine: 
                        stopOnMatch = True
                    else:
                        stopOnMatch = False
        if quitOnMatchRe is not None:
            if isinstance(stopOnMatchRe, str):
                pass
            else:
                if quitOnMatchRe.findall(lastLine):
                    stopOnMatch = True
        if tableEndRe is not None:
            if tableEndRe.findall(lastLine):
                stopOnMatch = True
        if stopOnMatch:
            return True, parserDict, globalDict, localDict
        else:
            # Check parser options to determine 
            # the table properties
            header = False 
            wrap = False 
            storeLines = 0
            lineFilter = None 
            headerList = None 
            headerListCheck = None 
            headersave = None 
            skipLines = 0
            if "header" in parserOptions.keys():
                header = parserOptions["header"] if (parserOptions["header"] is not None or 
                        parserOptions["header"] is not False) else False
            if "headersave" in parserOptions.keys():
                if parserOptions["headersave"]:
                    try:
                        headersave = getattr(self, parserOptions["headersave"]) 
                    except(AttributeError):
                        headersave = {}
                        setattr(self, parserOptions["headersave"], headersave)
                else:
                    headersave = None
            if "headerList" in parserOptions.keys():
                if parserOptions["headerList"]:
                    headerListCheck = getattr(self, parserOptions["headerList"]) 
                    if headerListCheck is not None:
                        headerList = headerListCheck
                    else:
                        headerList = parserOptions["headerList"]
                else:
                    headerList = None
            if "wrap" in parserOptions.keys():
                wrap = parserOptions["wrap"] if parserOptions["wrap"] else False
            if "tablelines" in parserOptions.keys():
                storeLines = parserOptions["tablelines"] if parserOptions["tablelines"]>0 else 0
            if "lineFilter" in parserOptions.keys():
                lineFilter = parserOptions["lineFilter"] if isinstance(parserOptions["lineFilter"], dict) else None
            if "skiplines" in parserOptions.keys():
                skipLines = parserOptions["skiplines"]
            if (skipLines>0 and skipLines>parserDict["firstLine"]):
                linenum = -1
            else:
                linenum = parserDict["numStoredLines"]
                if linenum < storeLines:
                    parserDict.update({"storedLines" : parserDict["storedLines"] + conv_str(lastLine)})
                    linenum += 1
                    parserDict.update({"numStoredLines" : linenum})
            if storeLines == 0:
                parserDict.update({"storedLines" : conv_str(lastLine)})
                linenum = 0
                parserDict.update({"numStoredLines" : linenum})

            mNameDict = getattr(self, matchNameDict)

            if linenum < storeLines:
                #linenum += 1
                #parserDict.update({"numStoredLines" : linenum})
                pass
            elif linenum == storeLines:
                storedLines = parserDict["storedLines"]
                storedText = ''
                # Convert with filter
                if lineFilter is not None:
                    for line in storedLines.splitlines():
                        for k, v in lineFilter.items():
                            line = line.replace(k, v)
                        storedText = storedText + line + '\n'
                else:
                    storedText =  storedLines
                # Extract table header and/or values
                lcount = 0
                hlist = []
                vlist = []
                for line in storedText.splitlines():
                    if header: 
                        if lcount<1:
                            hlist.extend(line.split())
                        else:
                            vlist.extend(line.split())
                    else:
                        if headerList is not None:
                            if isinstance(headerList, (tuple, list)):
                                for hitem in headerList:
                                    hlist.append(hitem)
                            elif isinstance(headerList, dict):
                                # Making sure the keys will be copied in ordered as it is assigned
                                for hi in range(0, len(headerList.keys())):
                                    for hitem, order in headerList.items():
                                        if hi == order:
                                            hlist.append(hitem)
                            elif isinstance(headerList, str):
                                for hitem in headerList.split():
                                    hlist.append(hitem)
                        vlist.extend(line.split())
                    lcount += 1
                # Reconvert to original
                if lineFilter is not None:
                    tmplist = []
                    for hitem in hlist:
                        for k, v in lineFilter.items():
                            hitem = hitem.replace(v, k)
                        tmplist.append(hitem)
                    hlist = tmplist
                #print("PRINTING table_parser lcount:hlist:vlist:",lcount,hlist,vlist)
                anythingmatched = False
                # Search for dictionary keys
                if vlist:
                    for cName, key in getDict_MetaStrInDict(mNameDict).items():
                        if onlyCaseSensitive is not True:
                            cName=cName.upper()
                            hlist=[h.upper() for h in hlist]
                        if cName in hlist:
                            #print("PRINTING table_parser save:",cName,v)
                            hindex = hlist.index(cName)
                            if hindex < len(vlist):
                                v=vlist[hlist.index(cName)]
                                if key in list(parser.lastMatch.keys()):
                                    parser.lastMatch[key]=v
                                else:
                                    mNameDict[key].value=v
                                    mNameDict[key].activeInfo=True
                                anythingmatched = True
                else:
                    if(header is not None and headersave is not None):
                        hlistc=0
                        # Saving keys with their orders in the list
                        for cName in hlist:
                            headersave.update({cName : int(hlistc)})
                            hlistc+=1
                        anythingmatched = True
                if anythingmatched:
                    #print("PRINTING table_parser header:",headersave)
                    # We have matched keywords to update at sections
                    # if the active sections are defined
                    cntrlDict = None
                    cntrlsec = None
                    cntrlsave = None
                    if "controldict" in parserOptions:
                        cntrlDict = getattr(self,parserOptions["controldict"])
                    if "controlsections" in parserOptions:
                        cntrlsec = parserOptions["controlsections"]
                    if "controlsave" in parserOptions:
                        if cntrlDict:
                            if parserOptions["controlsave"] in cntrlDict:
                                secDict = {}
                                cntrlsave = cntrlDict[parserOptions["controlsave"]]
                                if isinstance(cntrlsave,dict):
                                    secDict = cntrlsave
                                if cntrlsec is not None:
                                    for sec in cntrlsec:
                                        secDict.update({sec:True})
                                else:
                                    addname = "table_store"
                                    if parserDict["parserID"] is not None:
                                        addname = addname + str(parserDict["parserID"])
                                    secDict.update({addname:True})
                                cntrlDict.update({parserOptions["controlsave"]:secDict})
                        else:
                            secDict = {}
                            cntrlsave = getattr(self,parserOptions["controlsave"])
                            if isinstance(cntrlsave,dict):
                                secDict = cntrlsave
                            if cntrlsec is not None:
                                for sec in cntrlsec:
                                    secDict.update({sec:True})
                            else:
                                addname = "table_store"
                                if parserDict["parserID"] is not None:
                                    addname = addname + str(parserDict["parserID"])
                                secDict.update({addname:True})
                            setattr(self,parserOptions["controlsave"],secDict)
                        setattr(self,parserOptions["controldict"],cntrlDict)

                linenum = 0
                parserDict.update({"numStoredLines" : linenum})
                parserDict.update({"storedLines"    : ''})
                parserDict.update({"parserSuccess"  : True})
            else:
                linenum = 0
                parserDict.update({"numStoredLines" : linenum})
                parserDict.update({"storedLines"    : ''})

            if(mNameDict is not None and updateMatchDict):
                setattr(self, matchNameDict, mNameDict)
            if(headersave is not None and anythingmatched is True):
                setattr(self, parserOptions["headersave"], headersave)
            return False, parserDict, globalDict, localDict

    def table_parser(self, parser, stopOnMatchStr, quitOnMatchStr, metaNameStart, matchNameList, 
            matchNameDict, updateMatchDict, onlyCaseSensitive, stopOnFirstLine, parserOptions, 
            entryline=None, parserID=None, parsername=None, globalDict=None, localDict=None, parent=None):
        record = None
        replay = None
        replayCount = None
        if parent is None:
            parent = 0
        rank = parent + 1
        if globalDict is not None:
            record = globalDict["record"]
            replay = globalDict["replay"]
            replayCount = globalDict["replayCount"]
        if record is None:
            record = False
        if replayCount is None:
            replayCount = 0
        if entryline is None:
            lastLine = parser.fIn.fInLine
        else:
            lastLine = entryline
        #print("PRINTING table_parser:",lastLine)
        parserDict = {
                "firstLine"   : 0,
                "storedLines" : '',
                "numStoredLines" : 0,
                "parserID" : parserID,
                "parserSuccess" : False,
                "parent" : parent,
                "rank" : rank
                }
        updateLastLine = False
        parsercontinue = True
        parsercntrlattr = None
        parsercntrlin = None
        lookupdict = None
        lookupvals = None
        maxLines = None
        if "movetostopline" in parserOptions:
            updateLastLine = parserOptions["movetostopline"]
        if "lookupdict" in parserOptions:
            lookupdict = getattr(self,parserOptions["lookupdict"])
        if "parsercntrlattr" in parserOptions:
            parsercntrlattr = getattr(self,parserOptions["parsercntrlattr"])
        if "parsercntrlin" in parserOptions:
            if lookupdict:
                parsercntrlin = parserOptions["parsercntrlin"]
                if parsercntrlin in lookupdict:
                    lookupvals = lookupdict[parsercntrlin]
            else:
                if isinstance(parserOptions["parsercntrlin"], str):
                    parsercntrlin = getattr(self,parserOptions["parsercntrlin"])
                else:
                    parsercntrlin = parserOptions["parsercntrlin"]
        if "maxtablelines" in parserOptions:
            maxLines = parserOptions["maxtablelines"] if parserOptions[
                    "maxtablelines"]>0 else None

        if(parsercntrlattr is not None and
            parsercntrlin is not None):
            if lookupvals is not None:
                if parsercntrlattr in lookupvals:
                    parsercontinue = True
                else:
                    parsercontinue = False

        if parsercontinue:
            # Continue search and store until the line matches with stopOnMatch.
            stopOnMatchRe = re.compile(stopOnMatchStr)
            quitOnMatchRe = None
            if quitOnMatchStr is not None:
                quitOnMatchRe = re.compile(quitOnMatchStr)
            rtn, parserDict, globalDict, localDict = self.table_store(parser, lastLine, 
                    stopOnMatchRe, quitOnMatchRe,
                    metaNameStart, matchNameList, 
                    matchNameDict, updateMatchDict, 
                    onlyCaseSensitive, stopOnFirstLine, 
                    parserDict, parserOptions, globalDict, localDict) 
            if maxLines is not None:
                if maxLines >= parserDict["numStoredLines"]:
                    rtn = True
            if rtn is not True:
                while True:
                    if record and replayCount>0:
                        lastLine = self.peekrecord()
                    else:
                        lastLine = self.peekline(parser)
                    parserDict.update({"firstLine" : parserDict["firstLine"] + 1})
                    if not lastLine:
                        break
                    else:
                        # Matched with stopOnMatch. Discarding the line and 
                        # return SimpleMatcher context. Can this line be 
                        # discarded since it is the end of line for input 
                        # control variables or end of namelist ?
                        rtn, parserDict, globalDict, localDict = self.table_store(parser, lastLine, 
                                stopOnMatchRe, quitOnMatchRe, 
                                metaNameStart, matchNameList, 
                                matchNameDict, updateMatchDict, 
                                onlyCaseSensitive, stopOnFirstLine, 
                                parserDict, parserOptions, globalDict, localDict)
                        if maxLines is not None:
                            if maxLines >= parserDict["numStoredLines"]:
                                rtn = True
                        if rtn:
                            if updateLastLine:
                                # Playing from record?
                                if record and replayCount>0:
                                    lastLine = self.recordList.readline()
                                else:
                                    lastLine = parser.fIn.readline()
                                    # If not replaying, are we recording?
                                    if record:
                                        self.recordList.write(lastLine)
                            break
                        else:
                            # Playing from record?
                            if record and replayCount>0:
                                lastLine = self.recordList.readline()
                            else:
                                lastLine = parser.fIn.readline()
                                # If not replaying, are we recording?
                                if record:
                                    self.recordList.write(lastLine)
        return lastLine, parserDict["parserSuccess"], globalDict, localDict

    def namelist_store(self, parser, lastLine, stopOnMatchRe, quitOnMatchRe, 
            metaNameStart, matchNameList, matchNameDict, updateMatchDict, 
            onlyCaseSensitive, stopOnFirstLine, parserDict, parserOptions,
            globalDict=None, localDict=None):
        skipThis = False
        if "skipThis" in parserDict:
            if parserDict["skipThis"] is True:
                skipThis = True
        record = None
        replay = None
        replayCount = None
        if globalDict is not None:
            record = globalDict["record"]
            replay = globalDict["replay"]
            replayCount = globalDict["replayCount"]
        if record is None:
            record = False
        if replayCount is None:
            replayCount = 0
        stopOnMatch = False
        if stopOnMatchRe.findall(lastLine):
            stopOnMatch = True
            if self.firstLine==0:
                if stopOnFirstLine: 
                    stopOnMatch = True
                else:
                    stopOnMatch = False
        if quitOnMatchRe is not None:
            if quitOnMatchRe.findall(lastLine):
                stopOnMatch = True
        if stopOnMatch:
            if "skipThis" in parserDict:
                if parserDict["skipThis"] is True and skipThis is True:
                    skipThis = False
                    parserDict["skipThis"] = False
            return True, parserDict, globalDict, localDict
        else:
            anythingmatched = False
            if skipThis is False:
                # If there is at least one namelist in the line, 
                # search for all others in the dictionary.
                mNameDict = getattr(self, matchNameDict)
                for cName, key in getDict_MetaStrInDict(mNameDict).items():
                    if onlyCaseSensitive is not True:
                        c2Name = cName.upper()
                        c3Name = cName.lower()
                        cText = "\s*%s|^%s|,%s" % (cName, cName, cName)
                        cText = cText + "|\s*%s|^%s|,%s" % (c2Name, c2Name, c2Name)
                        cText = cText + "|\s*%s|^%s|,%s" % (c3Name, c3Name, c3Name)
                    else:
                        cText = "\s*%s|^%s|,%s" % (cName, cName, cName)
                    if mNameDict[key].matchFirst:
                        matchFirst = int(mNameDict[key].matchFirst)
                        for numletters in range(matchFirst,len(cName)):
                            n1Name = cName[0:numletters]
                            n2Name = n1Name.upper()
                            n3Name = n1Name.lower()
                            cText = cText + "|\s*%s|^%s|,%s" % (n1Name, n1Name, n1Name)
                            cText = cText + "|\s*%s|^%s|,%s" % (n2Name, n2Name, n2Name)
                            cText = cText + "|\s*%s|^%s|,%s" % (n3Name, n3Name, n3Name)
                    if mNameDict[key].appendMatchesList:
                        appendMatchesList = mNameDict[key].appendMatchesList
                    else:
                        appendMatchesList = None
                    if mNameDict[key].appendMatchesUntil:
                        appendMatchesUntil = mNameDict[key].appendMatchesUntil
                    else:
                        appendMatchesUntil = None
                    if mNameDict[key].appendLinesUntil:
                        appendLinesUntil = mNameDict[key].appendLinesUntil
                    else:
                        appendLinesUntil = None
                    if mNameDict[key].appendMaxLines:
                        appendMaxLines = mNameDict[key].appendMaxLines
                    else:
                        appendMaxLines = None
                    if mNameDict[key].matchAlso:
                        matchAlso = mNameDict[key].matchAlso
                    else:
                        matchAlso = None
                    if mNameDict[key].matchWith:
                        matchWith = mNameDict[key].matchWith
                    else:
                        matchWith = ''
                    if mNameDict[key].replaceDict:
                        replaceDict = mNameDict[key].replaceDict
                    else:
                        replaceDict = None
                    if mNameDict[key].subFunc:
                        if isinstance(mNameDict[key].subFunc, str):
                            subFunc = getattr(self,mNameDict[key].subFunc)
                        else:
                            subFunc = mNameDict[key].subFunc
                    else:
                        subFunc = None
                    if mNameDict[key].addAsList:
                        addAsList = mNameDict[key].addAsList
                    else:
                        addAsList = None
                    if mNameDict[key].appendToList:
                        appendToList = getattr(self,mNameDict[key].appendToList)
                    else:
                        appendToList = None
                    matchThisParsy = namelist_matcher(cText, key, matchWith)
                    reDict={key:value for value in matchThisParsy.findall(lastLine)}
                    if matchAlso:
                        if matchAlso not in lastLine:
                            reDict=None
                    if reDict:
                        for k,v in reDict.items():
                            if k == key: 
                                if isinstance(v, str):
                                    if replaceDict is not None:
                                        for repK, repV in replaceDict.items():
                                            v=v.replace(repK, repV)
                                if subFunc is not None:
                                    v=subFunc(v)
                                if k in list(parser.lastMatch.keys()):
                                    parser.lastMatch[k]=v
                                else:
                                    if addAsList is not None:
                                        if mNameDict[k].value is None:
                                            mNameDict[k].value=[v]
                                        else:
                                            if isinstance(mNameDict[k].value, list):
                                                mNameDict[k].value.append(v)
                                            else:
                                                mNameDict[k].value=v
                                    else:
                                        mNameDict[k].value=v
                                    mNameDict[k].activeInfo=True
                                parserDict.update({"parserSuccess"  : True})
                                anythingmatched = True
                        if appendToList is not None:
                            setattr(self, mNameDict[key].appendToList, appendToList)
            if anythingmatched:
                # We have matched keywords to update at sections
                # if the active sections are defined
                cntrlDict = None
                cntrlsec = None
                cntrlsave = None
                if "controldict" in parserOptions:
                    cntrlDict = getattr(self,parserOptions["controldict"])
                if "controlsections" in parserOptions:
                    cntrlsec = parserOptions["controlsections"]
                if "controlsave" in parserOptions:
                    if cntrlDict:
                        if parserOptions["controlsave"] in cntrlDict:
                            secDict = {}
                            cntrlsave = cntrlDict[parserOptions["controlsave"]]
                            if isinstance(cntrlsave,dict):
                                secDict = cntrlsave
                            if cntrlsec is not None:
                                for sec in cntrlsec:
                                    secDict.update({sec:True})
                            else:
                                addname = "namelist_store"
                                if parserDict["parserID"] is not None:
                                    addname = addname + str(parserDict["parserID"])
                                secDict.update({addname:True})
                            cntrlDict.update({parserOptions["controlsave"]:secDict})
                    else:
                        secDict = {}
                        cntrlsave = getattr(self,parserOptions["controlsave"])
                        if isinstance(cntrlsave,dict):
                            secDict = cntrlsave
                        if cntrlsec is not None:
                            for sec in cntrlsec:
                                secDict.update({sec:True})
                        else:
                            addname = "namelist_store"
                            if parserDict["parserID"] is not None:
                                addname = addname + str(parserDict["parserID"])
                            secDict.update({addname:True})
                        setattr(self,parserOptions["controlsave"],secDict)
                if "controldict" in parserOptions:
                    setattr(self,parserOptions["controldict"],cntrlDict)

            if(mNameDict is not None and updateMatchDict):
                setattr(self, matchNameDict, mNameDict)
            if "skipThis" in parserDict:
                if parserDict["skipThis"] is True and skipThis is True:
                    skipThis = False
                    parserDict["skipThis"] = False
            return False, parserDict, globalDict, localDict

    def namelist_parser(self, parser, stopOnMatchStr, quitOnMatchStr, metaNameStart, matchNameList, 
            matchNameDict, updateMatchDict, onlyCaseSensitive, stopOnFirstLine, parserOptions, 
            entryline=None, parserID=None, parsername=None, globalDict=None, localDict=None, parent=None):
        record = None
        replay = None
        replayCount = None
        if parent is None:
            parent = 0
        rank = parent + 1
        if globalDict is not None:
            record = globalDict["record"]
            replay = globalDict["replay"]
            replayCount = globalDict["replayCount"]
        if record is None:
            record = False
        if replayCount is None:
            replayCount = 0
        if entryline is None:
            lastLine = parser.fIn.fInLine
        else:
            lastLine = entryline
        #print("PRINTING: name_list:",lastLine)
        parserDict = {
                "firstLine"   : 0,
                "storedLines" : lastLine,
                "parserID" : parserID,
                "parserSuccess" : False,
                "skipThis" : False,
                "oneLining" : False,
                "oneLineId" : None,
                "oneLine" : None,
                "saveLastLine" : None,
                "parent" : parent,
                "rank" : rank
                }
        self.firstLine = 0
        updateLastLine = False
        resetMatchNameDict = False
        oneLinerList = None
        # Check the captured line has Fortran namelist variables and store them.
        # Continue search and store until the line matches with stopOnMatch.
        if "movetostopline" in parserOptions:
            updateLastLine = parserOptions["movetostopline"]
        if "resetNameDict" in parserOptions:
            resetMatchNameDict = parserOptions["resetNameDict"]
        if "makeOneLinerBetween" in parserOptions:
            oneLinerList = parserOptions["makeOneLinerBetween"]
        if resetMatchNameDict is True:
            mNameDict = getattr(self, matchNameDict)
            if mNameDict is not None:
                for key, value in mNameDict.items():
                    mNameDict[key].value=None
                    mNameDict[key].activeInfo=False
                setattr(self, matchNameDict, mNameDict)
        stopOnMatchRe = re.compile(stopOnMatchStr)
        quitOnMatchRe = None
        if quitOnMatchStr is not None:
            quitOnMatchRe = re.compile(quitOnMatchStr)
        parserDict["saveLastLine"]=lastLine
        if oneLinerList is not None:
            for mId, mItem in enumerate(oneLinerList):
                sRe=re.compile(mItem[0])
                if sRe.findall(lastLine):
                    #print("PRINTING: oneLinerId, lastLine:", parserDict["oneLineId"], lastLine)
                    if parserDict["oneLining"] is False:
                        parserDict["oneLining"]=True
                    parserDict["oneLineId"]=mId
                    #print("PRINTING: oneLinerId, oneLine:", parserDict["oneLineId"], parserDict["oneLine"])
                    break
        if parserDict["oneLining"] is True:
            mItem=oneLinerList[parserDict["oneLineId"]]
            repOn=False
            if len(mItem)>2:
                if isinstance(mItem[2], list):
                    if len(mItem[2])>0:
                        if len(mItem[2][0])>0:
                            repOn=True
            if parserDict["oneLine"] is None:
                if repOn is True:
                    repLastLine=lastLine
                    for repItem in mItem[2]:
                        if len(repItem)<2:
                            repStr=""
                        else:
                            repStr=repItem[1]
                        repLastLine=repLastLine.replace(repItem[0], repStr)
                    parserDict["oneLine"]=repLastLine
                else:
                    parserDict["oneLine"]=lastLine
            else:
                if repOn is True:
                    repLastLine=lastLine
                    for repItem in mItem[2]:
                        if len(repItem)<2:
                            repStr=""
                        else:
                            repStr=repItem[1]
                        repLastLine=repLastLine.replace(repItem[0], repStr)
                    parserDict["oneLine"]=parserDict["oneLine"]+repLastLine
                else:
                    parserDict["oneLine"]=parserDict["oneLine"]+lastLine
        if parserDict["oneLining"] is True:
            skipThis = True
            fRe=re.compile(oneLinerList[parserDict["oneLineId"]][1])
            if fRe.findall(lastLine):
                skipThis = False
                lastLine = parserDict["oneLine"]
                parserDict["oneLine"] = None
                parserDict["oneLining"] = False
                parserDict["oneLineId"] = None
        rtn, parserDict, globalDict, localDict = self.namelist_store(parser, lastLine, 
                stopOnMatchRe, quitOnMatchRe,
                metaNameStart, matchNameList, 
                matchNameDict, updateMatchDict, 
                onlyCaseSensitive, stopOnFirstLine, 
                parserDict, parserOptions, globalDict, localDict)
        lastLine = parserDict["saveLastLine"]
        if rtn is not True:
            while True:
                if record and replayCount>0:
                    lastLine = self.peekrecord()
                else:
                    lastLine = self.peekline(parser)
                self.firstLine += 1
                if not lastLine:
                    break
                else:
                    # Matched with stopOnMatch. Discarding the line and return SimpleMatcher context.
                    # Can this line be discarded since it is the end of line for input control
                    # variables or end of namelist ?
                    parserDict["saveLastLine"]=lastLine
                    if oneLinerList is not None:
                        for mId, mItem in enumerate(oneLinerList):
                            sRe=re.compile(mItem[0])
                            if sRe.findall(lastLine):
                                #print("PRINTING: oneLinerId, lastLine:", parserDict["oneLineId"], lastLine)
                                if parserDict["oneLining"] is False:
                                    parserDict["oneLining"]=True
                                parserDict["oneLineId"]=mId
                                break
                        if parserDict["oneLining"] is True:
                            mItem=oneLinerList[parserDict["oneLineId"]]
                            repOn=False
                            if len(mItem)>2:
                                if isinstance(mItem[2], list):
                                    if len(mItem[2])>0:
                                        if len(mItem[2][0])>0:
                                            repOn=True
                            if parserDict["oneLine"] is None:
                                if repOn is True:
                                    repLastLine=lastLine
                                    for repItem in mItem[2]:
                                        if len(repItem)<2:
                                            repStr=""
                                        else:
                                            repStr=repItem[1]
                                        repLastLine=repLastLine.replace(repItem[0], repStr)
                                    parserDict["oneLine"]=repLastLine
                                else:
                                    parserDict["oneLine"]=lastLine
                            else:
                                if repOn is True:
                                    repLastLine=lastLine
                                    for repItem in mItem[2]:
                                        if len(repItem)<2:
                                            repStr=""
                                        else:
                                            repStr=repItem[1]
                                        repLastLine=repLastLine.replace(repItem[0], repStr)
                                    parserDict["oneLine"]=parserDict["oneLine"]+repLastLine
                                else:
                                    parserDict["oneLine"]=parserDict["oneLine"]+lastLine
                    if parserDict["oneLining"] is True:
                        skipThis = True
                        fRe=re.compile(oneLinerList[parserDict["oneLineId"]][1])
                        if fRe.findall(lastLine):
                            skipThis = False
                            lastLine = parserDict["oneLine"]
                            #print("PRINTING: oneLinerId, oneLine:", parserDict["oneLineId"], parserDict["oneLine"])
                            parserDict["oneLine"] = None
                            parserDict["oneLining"] = False
                            parserDict["oneLineId"] = None
                    rtn, parserDict, globalDict, localDict = self.namelist_store(parser, lastLine, 
                            stopOnMatchRe, quitOnMatchRe, 
                            metaNameStart, matchNameList, 
                            matchNameDict, updateMatchDict, 
                            onlyCaseSensitive, stopOnFirstLine, 
                            parserDict, parserOptions, globalDict, localDict)
                    lastLine = parserDict["saveLastLine"]
                    if rtn:
                        if updateLastLine:
                            # Playing from record?
                            if record and replayCount>0:
                                lastLine = self.recordList.readline()
                            else:
                                lastLine = parser.fIn.readline()
                                # If not replaying, are we recording?
                                if record:
                                    self.recordList.write(lastLine)
                        break
                    else:
                        # Playing from record?
                        if record and replayCount>0:
                            lastLine = self.recordList.readline()
                        else:
                            lastLine = parser.fIn.readline()
                            # If not replaying, are we recording?
                            if record:
                                self.recordList.write(lastLine)
        return lastLine, parserDict["parserSuccess"], globalDict, localDict

    def check_namelist_store(self, parser, lastLine, stopOnMatchRe, quitOnMatchRe, metaNameStart, 
            matchNameList, matchNameDict, onlyCaseSensitive, stopOnFirstLine, parserDict):
        stopOnMatch = False
        if stopOnMatchRe.findall(lastLine):
            stopOnMatch = True
            if self.firstLine==0:
                if stopOnFirstLine: 
                    stopOnMatch = True
                else:
                    stopOnMatch = False
        if quitOnMatchRe is not None:
            if quitOnMatchRe.findall(lastLine):
                stopOnMatch = True
        if stopOnMatch:
            return True
        else:
            # If there is at least one namelist in the line, 
            # search for all others in the dictionary.
            if self.MD is not True:
                newLine = parser.fIn.readline()
                lastLine = ' = '.join([ "%s" % str(line) for line in zip(lastLine, newLine)])
            for cName, key in getDict_MetaStrInDict(matchNameDict).items():
                reDict={key:value for value in 
                        re.compile(r"(?:\s%s|^%s|,%s)\s*=\s*(?:'|\")?"
                                    "(?P<%s>[\-+0-9.a-zA-Z:]+)(?:'|\")?\s*,?" 
                        % (cName, cName, cName, key)).findall(lastLine)}
                if onlyCaseSensitive is not True:
                    reDict.update({key:value for value in 
                                   re.compile(r"(?:\s%s|^%s|,%s)\s*=\s*(?:'|\")?"
                                               "(?P<%s>[\-+0-9.a-zA-Z:]+)(?:'|\")?\s*,?" 
                                   % (cName.upper(), cName.upper(), 
                                      cName.upper(), key)).findall(lastLine)})
                if reDict:
                    for k,v in reDict.items():
                        if k == key: 
                            if k in list(parser.lastMatch.keys()):
                                parser.lastMatch[k]=v
                            else:
                                matchNameDict[k].value=v
                                matchNameDict[k].activeInfo=True
            return False

    def adHoc_read_namelist_stop_parsing(self, parser, stopOnMatchStr, quitOnMatchStr, 
            metaNameStart, matchNameList, matchNameDict, onlyCaseSensitive, stopOnFirstLine):
        lastLine = parser.fIn.fInLine
        parserDict = {
                "firstLine"   : 0,
                "storedLines" : lastLine,
                }
        self.firstLine = 0
        # Check the captured line has Fortran namelist variables and store them.
        # Continue search and store until the line matches with stopOnMatch.
        stopOnMatchRe = re.compile(stopOnMatchStr)
        quitOnMatchRe = None
        if quitOnMatchStr is not None:
            quitOnMatchRe = re.compile(quitOnMatchStr)
        if self.check_namelist_store(parser, lastLine, 
                stopOnMatchRe, quitOnMatchRe,
                metaNameStart, matchNameList, 
                matchNameDict, onlyCaseSensitive, 
                stopOnFirstLine, parserDict) is not True:
            while True:
                lastLine = self.peekline(parser)
                self.firstLine += 1
                if not lastLine:
                    break
                else:
                    # Matched with stopOnMatch. Discarding the line and return SimpleMatcher context.
                    # Can this line be discarded since it is the end of line for input control
                    # variables or end of namelist ?
                    if self.check_namelist_store(parser, lastLine, 
                            stopOnMatchRe, quitOnMatchRe, 
                            metaNameStart, matchNameList, 
                            matchNameDict, onlyCaseSensitive,
                            stopOnFirstLine, parserDict):
                        break
                    else:
                        lastLine = parser.fIn.readline()

    def check_commandline_args(self, parser, lastLine, stopOnMatchRe, quitOnMatchRe, 
            metaNameStart, matchNameList, matchNameDict, onlyCaseSensitive, stopOnFirstLine,
            parserDict, commandLineMatchRe, commandLineFunction):
        stopOnMatch = False
        if stopOnMatchRe.findall(lastLine):
            stopOnMatch = True
            if parserDict["firstLine"]==0:
                if stopOnFirstLine: 
                    stopOnMatch = True
                else:
                    stopOnMatch = False
        if quitOnMatchRe is not None:
            if quitOnMatchRe.findall(lastLine):
                stopOnMatch = True
        # If there is at least one namelist in the line, 
        # search for all others in the dictionary.
        checkDict = getDict_MetaStrInDict(matchNameDict, matchNameList)
        cmdMatch = commandLineMatchRe.findall(lastLine)
        if cmdMatch:
            cmdDict = commandLineFunction(cmdMatch[0])
            for cName, key in checkDict.items():
                cmdName = cName.replace("-", "")
                if cmdName in cmdDict:
                    for k,v in cmdDict.items():
                        if k == cmdName: 
                            if k in list(parser.lastMatch.keys()):
                                parser.lastMatch[key]=v
                            else:
                                matchNameDict[key].value=v
                                matchNameDict[key].activeInfo=True
        if stopOnMatch:
            return True
        else:
            return False

    def adHoc_read_commandline_stop_parsing(self, parser, stopOnMatchStr, quitOnMatchStr, 
            metaNameStart, matchNameList, matchNameDict, onlyCaseSensitive, stopOnFirstLine,
            commandLineMatchStr, commandLineFunction):
        lastLine = parser.fIn.fInLine
        parserDict = {
                "firstLine"   : 0,
                "storedLines" : lastLine,
                }
        parserDict["firstLine"] = 0
        # Check the captured line has Fortran namelist variables and store them.
        # Continue search and store until the line matches with stopOnMatch.
        stopOnMatchRe = re.compile(stopOnMatchStr)
        commandLineMatchRe = re.compile(commandLineMatchStr)
        quitOnMatchRe = None
        if quitOnMatchStr is not None:
            quitOnMatchRe = re.compile(quitOnMatchStr)
        if self.check_commandline_args(parser, lastLine, 
                stopOnMatchRe, quitOnMatchRe,
                metaNameStart, matchNameList, 
                matchNameDict, onlyCaseSensitive, 
                stopOnFirstLine, parserDict, 
                commandLineMatchRe, commandLineFunction) is not True:
            while True:
                lastLine = self.peekline(parser)
                parserDict.update({"firstLine" : parserDict["firstLine"] + 1})
                if not lastLine:
                    break
                else:
                    # Matched with stopOnMatch. Discarding the line and return SimpleMatcher context.
                    # Can this line be discarded since it is the end of line for input control
                    # variables or end of namelist ?
                    if self.check_commandline_args(parser, lastLine, 
                            stopOnMatchRe, quitOnMatchRe, 
                            metaNameStart, matchNameList, 
                            matchNameDict, onlyCaseSensitive,
                            stopOnFirstLine, parserDict, 
                            commandLineMatchRe, commandLineFunction):
                        break
                    else:
                        lastLine = parser.fIn.readline()

    def adHoc_find_textfile_store(self, parser, 
            commandLineMatchStr, commandLineFunction):
        lastLine = parser.fIn.fInLine
        parserDict = {
                "firstLine"   : 0,
                "storedLines" : lastLine,
                }
        parserDict["firstLine"] = 0
        # Check the captured line has Fortran namelist variables and store them.
        # Continue search and store until the line matches with stopOnMatch.
        stopOnMatchRe = re.compile(stopOnMatchStr)
        commandLineMatchRe = re.compile(commandLineMatchStr)
        quitOnMatchRe = None
        if quitOnMatchStr is not None:
            quitOnMatchRe = re.compile(quitOnMatchStr)
        if self.check_commandline_args(parser, lastLine, 
                stopOnMatchRe, quitOnMatchRe,
                metaNameStart, matchNameList, 
                matchNameDict, onlyCaseSensitive, 
                stopOnFirstLine, parserDict, 
                commandLineMatchRe, commandLineFunction) is not True:
            while True:
                lastLine = self.peekline(parser)
                parserDict.update({"firstLine" : parserDict["firstLine"] + 1})
                if not lastLine:
                    break
                else:
                    # Matched with stopOnMatch. Discarding the line and return SimpleMatcher context.
                    # Can this line be discarded since it is the end of line for input control
                    # variables or end of namelist ?
                    if self.check_commandline_args(parser, lastLine, 
                            stopOnMatchRe, quitOnMatchRe, 
                            metaNameStart, matchNameList, 
                            matchNameDict, onlyCaseSensitive,
                            stopOnFirstLine, parserDict, 
                            commandLineMatchRe, commandLineFunction):
                        break
                    else:
                        lastLine = parser.fIn.readline()

    def check_text_store(self, parser, lastLine, stopOnMatchRe, quitOnMatchRe, 
            metaNameStart, metaNameStore, matchNameList, matchNameDict, onlyCaseSensitive, 
            stopOnFirstLine, storeFirstLine, storeStopQuitLine, parserDict):
        stopOnMatch = False
        key = metaNameStore
        parserDict.update({"storedLines" : parserDict["storedLines"] + conv_str(lastLine)})
        if stopOnMatchRe.findall(lastLine):
            stopOnMatch = True
            if parserDict["firstLine"]==0:
                if stopOnFirstLine: 
                    stopOnMatch = True
                else:
                    stopOnMatch = False
        if quitOnMatchRe is not None:
            if quitOnMatchRe.findall(lastLine):
                stopOnMatch = True
        if stopOnMatch:
            if storeStopQuitLine:
                if key in list(parser.lastMatch.keys()):
                    parser.lastMatch[key]=parserDict["storedLines"]
                else:
                    parser.lastMatch.update({key : parserDict["storedLines"]})
                if matchNameDict:
                    if key in matchNameDict:
                        matchNameDict[key].value=parserDict["storedLines"]
                        matchNameDict[key].activeInfo=True
            return True
        else:
            # If parsing is not stop or quit,
            # store the parsed string in text
            # with the given meta info name.
            # This will also copy the value to 
            # matchNameDict if the metaNameStore 
            # exists in the dictionary keys.
            if(parserDict["firstLine"]==0 and storeFirstLine):
                parserDict.update({"storedLines" : conv_str(lastLine)})
                if key in list(parser.lastMatch.keys()):
                    parser.lastMatch[key]=parserDict["storedLines"]
                else:
                    parser.lastMatch.update({key : parserDict["storedLines"]})
                if matchNameDict:
                    if key in matchNameDict:
                        matchNameDict[key].value=parserDict["storedLines"]
                        matchNameDict[key].activeInfo=True
            else:
                if key in list(parser.lastMatch.keys()):
                    parser.lastMatch[key]=parserDict["storedLines"]
                else:
                    parser.lastMatch.update({key : parserDict["storedLines"]})
                if matchNameDict:
                    if key in matchNameDict:
                        matchNameDict[key].value=parserDict["storedLines"]
                        matchNameDict[key].activeInfo=True
            return False

    def adHoc_read_store_text_stop_parsing(self, parser, stopOnMatchStr, quitOnMatchStr, 
            metaNameStart, metaNameStore, matchNameList, matchNameDict, onlyCaseSensitive, 
            stopOnFirstLine, storeFirstLine, storeStopQuitLine, onQuitRunFunction):
        lastLine = parser.fIn.fInLine
        parserDict = {
                "firstLine"   : 0,
                "storedLines" : lastLine,
                }
        parserDict["firstLine"] = 0
        # Store all parsed lines until stopOnMatchStr or quitOnMatchStr is encountered.
        stopOnMatchRe = re.compile(stopOnMatchStr)
        quitOnMatchRe = None
        if quitOnMatchStr is not None:
            quitOnMatchRe = re.compile(quitOnMatchStr)
        # Check first line to store and stop/quit
        if self.check_text_store(parser, lastLine, 
                stopOnMatchRe, quitOnMatchRe,
                metaNameStart, metaNameStore, 
                matchNameList, matchNameDict, 
                onlyCaseSensitive, stopOnFirstLine, 
                storeFirstLine, storeStopQuitLine,
                parserDict) is not True:
            # Check all other lines to store and stop/quit
            while True:
                lastLine = self.peekline(parser)
                parserDict.update({"firstLine" : parserDict["firstLine"] + 1})
                if not lastLine:
                    break
                else:
                    # Matched with stopOnMatch. The line will be procesed and 
                    # parsing will continue by returning SimpleMatcher context.
                    if self.check_text_store(parser, lastLine, 
                            stopOnMatchRe, quitOnMatchRe, 
                            metaNameStart, metaNameStore, 
                            matchNameList, matchNameDict, 
                            onlyCaseSensitive, stopOnFirstLine, 
                            storeFirstLine, storeStopQuitLine,
                            parserDict):
                        break
                    else:
                        lastLine = parser.fIn.readline()
        if onQuitRunFunction is not None:
            onQuitRunFunction(parser)

    def check_subparsers(self, parser, lastLine, stopOnMatchRe, quitOnMatchRe, subParsers, 
            ordered, parserDict, secDict, onStartRunFunction, onQuitRunFunction, 
            onlySubParsersReadLine, maxcycle=None, parent=None):
        record = None
        replay = None
        replayCount = None
        if parserDict is not None:
            record = parserDict["record"]
            replay = parserDict["replay"]
            replayCount = parserDict["replayCount"]
        if record is None:
            record = False
        if replayCount is None:
            replayCount = 0
        if maxcycle is None:
            maxcycle=1000000
        cyclenum = 0
        if parent is None:
            parent = 0
        rank = parent
        stopOnEOF = True
        waitFirstCycle = False
        if "waitFirstCycle" in parserDict:
            waitFirstCycle = parserDict["waitFirstCycle"]
        if "peekedline" in secDict:
            secDict.update({"peekedline":None})
        while True:
            stopOnMatch = False
            if waitFirstCycle is not True:
                if stopOnMatchRe.findall(lastLine):
                    stopOnMatch = True
                if quitOnMatchRe is not None:
                    if quitOnMatchRe.findall(lastLine):
                        stopOnMatch = True
            if subParsers:
                if ordered:
                    for parser_id in range(len(subParsers)):
                        ismatched=False
                        if(subParsers[parser_id]["startReStr"] is None or 
                            "AlwaysMatch" in subParsers[parser_id]["startReStr"] or 
                            "AlwaysRun" in subParsers[parser_id]["startReStr"]):
                            ismatched=True
                        else:
                            currentMatchRe = re.compile(subParsers[parser_id]["startReStr"])
                            if currentMatchRe.findall(lastLine):
                                ismatched=True
                        if ismatched:
                            skipme = False
                            if "waitlist" in subParsers[parser_id]:
                                if subParsers[parser_id]["waitlist"]:
                                    for waitlist in subParsers[parser_id]["waitlist"]:
                                        waitcheck = [False for wl in waitlist]
                                        for wi in range(len(waitlist)):
                                            if waitlist[wi] in parserDict["parserNameList"]:
                                                waitcheck[wi] = parserDict["parserNameList"][waitlist[wi]]
                                        if False not in waitcheck:
                                            skipme = True
                                            break
                            if(parserDict["parserCheckList"][parser_id] is False and 
                                False not in parserDict["parserCheckList"][0:parser_id] and 
                                skipme is False):
                                parserSuccess = False
                                parserName = None
                                if "parsername" in subParsers[parser_id]:
                                    parserName=subParsers[parser_id]["parsername"]
                                currentParser = getattr(self, subParsers[parser_id]["parser"])
                                parserDict.update({"matchStr":subParsers[parser_id]["startReStr"]})
                                parserDict.update({"matchLine":lastLine})
                                lastLine, parserSuccess, parserDict, secDict = currentParser(parser,
                                        subParsers[parser_id]["stopOnMatchStr"],
                                        subParsers[parser_id]["quitOnMatchStr"],
                                        subParsers[parser_id]["metaNameStart"],
                                        subParsers[parser_id]["matchNameList"],
                                        subParsers[parser_id]["matchNameDict"],
                                        subParsers[parser_id]["updateMatchDict"],
                                        subParsers[parser_id]["onlyCaseSensitive"],
                                        subParsers[parser_id]["stopOnFirstLine"],
                                        subParsers[parser_id]["parserOptions"],
                                        entryline=lastLine, parserID=parser_id,
                                        parsername=parserName,
                                        globalDict=parserDict,
                                        localDict=secDict,
                                        parent=rank
                                        )
                                parserDict["parserCheckList"][parser_id] = True
                                if "parsername" in subParsers[parser_id]:
                                    parserDict["parserNameList"].update({
                                        subParsers[parser_id]["parsername"] : parserSuccess
                                        })
                        if False not in parserDict["parserCheckList"]:
                            # all parsers are matched, reseting all
                            parserDict.update({"parserCheckList" : [False for i in range(len(subParsers))]})
                else:
                    parser_num = 0
                    for subParser in subParsers:
                        ismatched=False
                        if(subParser["startReStr"] is None or 
                            "AlwaysMatch" in subParser["startReStr"] or
                            "AlwaysRun" in subParser["startReStr"]):
                            ismatched=True
                        else:
                            currentMatchRe = re.compile(subParser["startReStr"])
                            if currentMatchRe.findall(lastLine):
                                ismatched=True
                        if ismatched:
                            skipme = False
                            if "waitlist" in subParser:
                                if subParser["waitlist"]:
                                    for waitlist in subParser["waitlist"]:
                                        waitcheck = [False for wl in waitlist]
                                        for wi in range(len(waitlist)):
                                            if waitlist[wi] in parserDict["parserNameList"]:
                                                waitcheck[wi] = parserDict["parserNameList"][waitlist[wi]]
                                        skipme = True if False in waitcheck else False
                                        if skipme is False:
                                            break
                            if skipme is False:
                                currentParser = getattr(self, subParser["parser"])
                                parserSuccess = False
                                parserName = None
                                if "parsername" in subParser:
                                    parserName=subParser["parsername"]
                                parserDict.update({"matchStr":subParser["startReStr"]})
                                parserDict.update({"matchLine":lastLine})
                                lastLine, parserSuccess, parserDict, secDict = currentParser(parser,
                                        subParser["stopOnMatchStr"],
                                        subParser["quitOnMatchStr"],
                                        subParser["metaNameStart"],
                                        subParser["matchNameList"],
                                        subParser["matchNameDict"],
                                        subParser["updateMatchDict"],
                                        subParser["onlyCaseSensitive"],
                                        subParser["stopOnFirstLine"],
                                        subParser["parserOptions"],
                                        entryline=lastLine, parserID=parser_num,
                                        parsername=parserName,
                                        globalDict=parserDict,
                                        localDict=secDict,
                                        parent=rank
                                        )
                                if "parsername" in subParser:
                                    parserDict["parserNameList"].update({
                                        subParser["parsername"] : parserSuccess
                                        })
                        parser_num += 1

            cyclenum += 1
            if waitFirstCycle is True:
                if stopOnMatchRe.findall(lastLine):
                    #print("PRINTING: check_subparser Rank waitFirstCycle StopMatched:",rank,stopOnMatchRe.findall(lastLine),lastLine)
                    stopOnMatch = True
                if quitOnMatchRe is not None:
                    if quitOnMatchRe.findall(lastLine):
                        #print("PRINTING: check_subparser Rank QuitMatched:",quitOnMatchRe.findall(lastLine),lastLine)
                        stopOnMatch = True
            if(onlySubParsersReadLine is None or
               onlySubParsersReadLine is False):
                break
            else:
                # if stop match is reached quit loop
                if stopOnMatch:
                    break
                # or if maximum loop cycle is 
                # reached than quit loop
                if cyclenum > maxcycle:
                    break
            if stopOnEOF:
                eofLine = None
                if record and replayCount>0:
                    eofLine = self.peekrecord()
                else:
                    eofLine = self.peekline(parser)
                if not eofLine:
                    stopOnMatch = True
                    break

        if stopOnMatch:
            return True, lastLine, parserDict, secDict
        else:
            return False, lastLine, parserDict, secDict

    def subparser_caller(self, parser, stopOnMatchStr, quitOnMatchStr, metaNameStart, 
            matchNameList, matchNameDict, updateMatchDict, onlyCaseSensitive, stopOnFirstLine, 
            parserOptions, entryline=None, parserID=None, parsername=None, globalDict=None,
            localDict=None, parent=None):
        record = None
        replay = None
        replayCount = None
        if parent is None:
            parent = 0
        rank = parent + 1
        if globalDict is not None:
            record = globalDict["record"]
            replay = globalDict["replay"]
            replayCount = globalDict["replayCount"]
        if record is None:
            record = False
        if replayCount is None:
            replayCount = 0
        if entryline is None:
            lastLine = parser.fIn.fInLine
        else:
            lastLine = entryline
        #print("PRINTING: subparser_caller,rank:",rank,lastLine)
        parserDict = {
                "firstLine"   : 0,
                "storedLines" : '',
                "numStoredLines" : 0,
                "parserID" : parserID,
                "parent" : parent,
                "rank" : rank
                }
        #parserDict.update({"firstLine" : parserDict["firstLine"] + 1})
        stopOnMatchRe = None
        quitOnMatchRe = None
        if stopOnMatchStr is not None:
            stopOnMatchRe = re.compile(stopOnMatchStr)
        if quitOnMatchStr is not None:
            quitOnMatchRe = re.compile(quitOnMatchStr)
        parserSuccess = False
        ordered=False
        onStartRunFunction=None
        onQuitRunFunction=None 
        onlySubParsersReadLine=False
        lookupdict = None
        lookupvals = None
        subparserMapper = None
        subparserSection = None
        sectionControlAttr = None
        matchFirstLetters = None
        sectionOnDict = {}
        if globalDict is not None:
            if "sectionOnDict" in globalDict:
                sectionOnDict = globalDict["sectionOnDict"]
        if "onStartRunFunction" in parserOptions:
            onStartRunFunction = parserOptions["onStartRunFunction"]
        if "onQuitRunFunction" in parserOptions:
            onQuitRunFunction = parserOptions["onQuitRunFunction"]
        if "ordered" in parserOptions:
            ordered = parserOptions["ordered"]
        if "onlySubParsersReadLine" in parserOptions:
            onlySubParsersReadLine = parserOptions["onlySubParsersReadLine"]
        if "subparserMapper" in parserOptions:
            subparserMapper = parserOptions["subparserMapper"]
        if "subparserSection" in parserOptions:
            subparserSection = parserOptions["subparserSection"]
        if "matchFirstLetters" in parserOptions:
            matchFirstLetters = parserOptions["matchFirstLetters"]
        if "waitFirstCycle" in parserOptions:
            if globalDict is not None:
                globalDict.update({"waitFirstCycle":parserOptions["waitFirstCycle"]})

        if onStartRunFunction is not None:
            onStartRunFunction(parser)

        sO = open_section
        #supB = parser.backend.superBackend
        supB = parser.backend
        sectionAct = None
        if subparserMapper is not None:
            for matchStr, listname in subparserMapper.items():
                cText = matchStr
                if onlyCaseSensitive is not True:
                    c2Name = matchStr.upper()
                    c3Name = matchStr.lower()
                    cText = "%s|%s|%s" % (matchStr, c2Name, c3Name)
                if matchFirstLetters is not None:
                    matchFirst = int(matchFirstLetters)
                    for numletters in range(matchFirst,len(matchStr)):
                        n1Name = cName[0:numletters]
                        n2Name = n1Name.upper()
                        n3Name = n1Name.lower()
                        cText = cText + "|%s|%s|%s" % (n1Name, n2Name, n3Name)
                matchRe = re.compile(cText)
                if matchRe.findall(lastLine):
                    subParsers = getattr(self,listname)
                    parserSuccess = True
                    if subparserSection is not None:
                        if matchStr in subparserSection:
                            sectionAct = subparserSection[matchStr]
                    if sectionAct is not None:
                        sectionOpenClose = sectionAct[1]
                        if not isinstance(sectionOpenClose, str):
                            sectionOpenClose = ''
                        sectionOn = getattr(self,sectionAct[2])
                        sectionGIndex = getattr(self,sectionAct[3])
                        if sectionAct[0] not in sectionOnDict:
                            sectionMatchCount = 0
                            sectionOnDict.update({sectionAct[0]:[sectionOn,sectionMatchCount]})
                        else:
                            sectionMatchCount = sectionOnDict[sectionAct[0]][1] + 1
                            sectionOnDict.update({
                                sectionAct[0]:[sectionOnDict[sectionAct[0]][0],sectionMatchCount]})
                        #print("PRINTING: subparser_caller, section Count:",sectionMatchCount)
                        #print("PRINTING: subparser_caller, section On:",sectionOn)
                        #print("PRINTING: subparser_caller, section First On:",sectionOnDict[sectionAct[0]][0])
                        if sectionOn is False:
                            if "AUTO" in sectionOpenClose.upper():
                                if sectionOnDict[sectionAct[0]][0] is False:
                                    with sO(supB, sectionAct[0]):
                                        rtn, lastLine, globalDict, parserDict = self.check_subparsers(
                                                parser, lastLine, stopOnMatchRe, quitOnMatchRe,
                                                subParsers, ordered, globalDict, parserDict, 
                                                onStartRunFunction, onQuitRunFunction, 
                                                onlySubParsersReadLine, parent=rank) 
                                else:
                                    sectionGIndex = supB.openSection(sectionAct[0])
                                    setattr(self,sectionAct[3], sectionGIndex)
                                    rtn, lastLine, globalDict, parserDict = self.check_subparsers(
                                            parser, lastLine, stopOnMatchRe, quitOnMatchRe,
                                            subParsers, ordered, globalDict, parserDict, 
                                            onStartRunFunction, onQuitRunFunction, 
                                            onlySubParsersReadLine, parent=rank) 
                            elif "OPEN" in sectionOpenClose.upper():
                                sectionGIndex = supB.openSection(sectionAct[0])
                                setattr(self,sectionAct[3], sectionGIndex)
                                rtn, lastLine, globalDict, parserDict = self.check_subparsers(
                                        parser, lastLine, stopOnMatchRe, quitOnMatchRe,
                                        subParsers, ordered, globalDict, parserDict, 
                                        onStartRunFunction, onQuitRunFunction, 
                                        onlySubParsersReadLine, parent=rank) 
                            else:
                                rtn, lastLine, globalDict, parserDict = self.check_subparsers(
                                        parser, lastLine, stopOnMatchRe, quitOnMatchRe,
                                        subParsers, ordered, globalDict, parserDict, 
                                        onStartRunFunction, onQuitRunFunction, 
                                        onlySubParsersReadLine, parent=rank) 
                        else:
                            if "AUTO" in sectionOpenClose.upper():
                                if sectionMatchCount>0:
                                    if sectionOnDict[sectionAct[0]][0] is True:
                                        supB.closeSection(sectionAct[0], sectionGIndex)
                                        sectionGIndex = supB.openSection(sectionAct[0])
                                        setattr(self,sectionAct[3], sectionGIndex)
                            rtn, lastLine, globalDict, parserDict = self.check_subparsers(
                                    parser, lastLine, stopOnMatchRe, quitOnMatchRe,
                                    subParsers, ordered, globalDict, parserDict, 
                                    onStartRunFunction, onQuitRunFunction, 
                                    onlySubParsersReadLine, parent=rank) 
                            if "CLOSE" in sectionOpenClose.upper():
                                supB.closeSection(sectionAct[0], sectionGIndex)
                            if "AUTO" in sectionOpenClose.upper():
                                if sectionOnDict[sectionAct[0]][0] is False:
                                    supB.closeSection(sectionAct[0], sectionGIndex)

                        # Infuse the correct gIndex to SimpleMatcher parsers to prevent 
                        # them fail before leaving the SmartParser.
                        #     Check if GIndex changed so far for the selected section
                        #     and update it if there is an open section left.
                        for SMid, SMcontext in enumerate(self.parser.context):
                            if sectionAct[0] in SMcontext.sections:
                                self.parser.context[SMid].sections[sectionAct[0]]=sectionGIndex
                    else:
                        rtn, lastLine, globalDict, parserDict = self.check_subparsers(
                                parser, lastLine, stopOnMatchRe, quitOnMatchRe,
                                subParsers, ordered, globalDict, parserDict, 
                                onStartRunFunction, onQuitRunFunction, 
                                onlySubParsersReadLine, parent=rank) 
                        parserDict.update({"firstLine" : parserDict["firstLine"] + 1})
                        if rtn:
                            break
                        else:
                            if onlySubParsersReadLine is False:
                                if record and replayCount>0:
                                    lastLine = self.recordList.readline()
                                else:
                                    lastLine = parser.fIn.readline()
                                    if record:
                                        self.recordList.write(lastLine)

#        stopOnMatch = False
#        if stopOnMatchRe.findall(lastLine):
#            stopOnMatch = True
#            if parserDict["firstLine"] == 0:
#                if stopOnFirstLine: 
#                    stopOnMatch = True
#                else:
#                    stopOnMatch = False
#        if quitOnMatchRe is not None:
#            if quitOnMatchRe.findall(lastLine):
#                stopOnMatch = True 

        if globalDict is not None:
            globalDict.update({"sectionOnDict":sectionOnDict})
       

        if onQuitRunFunction is not None:
            onQuitRunFunction(parser)

        return lastLine, parserSuccess, globalDict, localDict

    def adHoc_takingover_parsing(self, parser, stopOnMatchStr, quitOnMatchStr, stopControl=None, 
            record=False, replay=None, parseOnlyRecorded=False, subParsers=None, ordered=None, 
            onStartRunFunction=None, onQuitRunFunction=None, onlySubParsersReadLine=False, 
            onStartReplayRunFunction=None, onQuitReplayRunFunction=None, parent=None):
        if onStartRunFunction is not None:
            onStartRunFunction(parser)
        if ordered is None:
            ordered = False
        if stopControl is None:
            stopCntrl = True
        else:
            try:
                stopCntrl = getattr(self,stopControl)
            except(TypeError,ValueError,AttributeError):
                stopCntrl = True
        if record is not None or record is not False:
            record = True
        if replay is not None:
            try:
                replay = int(replay)
            except(ValueError,TypeError):
                replay = None
        lastLine = parser.fIn.fInLine
        if parent is None:
            parent = 0
        rank = parent
        parserDict = {
                "firstLine"  : 0,
                "parserCheckList" : [False for i in range(len(subParsers))],
                "parserNameList"  : {},
                "stopControl"     : stopCntrl,
                "record"          : record,
                "replay"          : replay,
                "replayCount"     : 0,
                "parent"          : parent
                }
        secDict = {}
        secDict.update(parserDict)
        # If true, we are recording the output in
        # every accurance of parser.fIn.readline
        #if record:
        #    self.recordList.write(lastLine)
        # parse lines until stopOnMatchStr or quitOnMatchStr is encountered.
        if onStartRunFunction is not None:
            onStartRunFunction(parser)
        stopOnMatchRe = None
        if stopOnMatchStr is not None:
            stopOnMatchRe = re.compile(stopOnMatchStr)
        quitOnMatchRe = None
        if quitOnMatchStr is not None:
            quitOnMatchRe = re.compile(quitOnMatchStr)
        if parseOnlyRecorded is False:
            # Check first line to store and stop/quit
            rtn, lastLine, parserDict, secDict = self.check_subparsers(
                parser, lastLine, stopOnMatchRe, quitOnMatchRe,
                subParsers, ordered, parserDict, secDict, 
                onStartRunFunction, onQuitRunFunction, 
                onlySubParsersReadLine, parent=rank) 
        else:
            rtn = False
            if stopOnMatchRe is not None:
                if stopOnMatchRe.findall(lastLine):
                    rtn = True
            if quitOnMatchRe is not None:
                if quitOnMatchRe.findall(lastLine):
                    rtn = True
        if secDict is not None:
            parserDict.update(secDict)
        if rtn is not True:
            # Check all other lines to store and stop/quit
            while True:
                lastLine = self.peekline(parser)
                parserDict.update({"firstLine" : parserDict["firstLine"] + 1})
                if not lastLine:
                    break
                else:
                    if parseOnlyRecorded is False:
                        # Matched with stopOnMatch. The line will be procesed and 
                        # parsing will continue by returning SimpleMatcher context.
                        rtn, lastLine, parserDict, secDict = self.check_subparsers(
                            parser, lastLine, stopOnMatchRe, quitOnMatchRe, 
                            subParsers, ordered, parserDict, secDict, 
                            onStartRunFunction, onQuitRunFunction, 
                            onlySubParsersReadLine, parent=rank)
                    else:
                        rtn = False
                        if stopOnMatchRe is not None:
                            if stopOnMatchRe.findall(lastLine):
                                rtn = True
                        if quitOnMatchRe is not None:
                            if quitOnMatchRe.findall(lastLine):
                                rtn = True
                    if secDict is not None:
                        parserDict.update(secDict)
                    if rtn:
                        if record:
                            self.recordList.write(lastLine)
                        break
                    else:
                        if onlySubParsersReadLine is False:
                            lastLine = parser.fIn.readline()
                            if record:
                                self.recordList.write(lastLine)
                        else:
                            if parseOnlyRecorded is True:
                                lastLine = parser.fIn.readline()
                                if record:
                                    self.recordList.write(lastLine)
        if record is True and replay is not None:
            if replay < 0:
                replayCount = 1
                while stopCntrl is not True:
                    if onStartReplayRunFunction is not None:
                        onStartReplayRunFunction[replayCount](parser)
                    self.recordList.seek(0)
                    lastLine = self.peekrecord()
                    parserDict["replayCount"]=replayCount
                    rtn, lastLine, parserDict, secDict = self.check_subparsers(
                            parser, lastLine, stopOnMatchRe, quitOnMatchRe,
                            subParsers, ordered, parserDict, secDict, 
                            onStartRunFunction, onQuitRunFunction, 
                            onlySubParsersReadLine, parent=rank) 
                    if secDict is not None:
                        parserDict.update(secDict)
                    if rtn is not True:
                        while True:
                            lastLine = self.peekrecord()
                            parserDict.update({"firstLine" : parserDict["firstLine"] + 1})
                            if not lastLine:
                                break
                            else:
                                rtn, lastLine, parserDict, secDict = self.check_subparsers(
                                        parser, lastLine, stopOnMatchRe, quitOnMatchRe, 
                                        subParsers, ordered, parserDict, secDict, 
                                        onStartRunFunction, onQuitRunFunction, 
                                        onlySubParsersReadLine, parent=rank)
                                if secDict is not None:
                                    parserDict.update(secDict)
                                if rtn:
                                    break
                                else:
                                    if onlySubParsersReadLine is False:
                                        lastLine = self.recordList.readline()
                    if onQuitReplayRunFunction is not None:
                        onQuitReplayRunFunction[replayCount](parser)
                    replayCount+=1
                    try:
                        stopCntrl = getattr(self,stopControl)
                    except(TypeError,ValueError,AttributeError):
                        stopCntrl = True
            elif replay > 0:
                replayCount = 1
                for repeat in range(replay):
                    # dirty function caller
                    if onStartReplayRunFunction is not None:
                        if replayCount in onStartReplayRunFunction:
                            if hasattr(self,onStartReplayRunFunction[replayCount]):
                                func = getattr(self,onStartReplayRunFunction[replayCount])
                                outfunc = func(parser)
                    self.recordList.seek(0)
                    lastLine = self.peekrecord()
                    parserDict["replayCount"]=repeat+1
                    rtn, lastLine, parserDict, secDict = self.check_subparsers(
                            parser, lastLine, stopOnMatchRe, quitOnMatchRe,
                            subParsers, ordered, parserDict, secDict, 
                            onStartRunFunction, onQuitRunFunction, 
                            onlySubParsersReadLine, parent=rank) 
                    if secDict is not None:
                        parserDict.update(secDict)
                    if rtn is not True:
                        while True:
                            lastLine = self.peekrecord()
                            parserDict.update({"firstLine" : parserDict["firstLine"] + 1})
                            if not lastLine:
                                break
                            else:
                                rtn, lastLine, parserDict, secDict = self.check_subparsers(
                                        parser, lastLine, stopOnMatchRe, quitOnMatchRe, 
                                        subParsers, ordered, parserDict, secDict, 
                                        onStartRunFunction, onQuitRunFunction, 
                                        onlySubParsersReadLine, parent=rank)
                                if secDict is not None:
                                    parserDict.update(secDict)
                                if rtn:
                                    break
                                else:
                                    if onlySubParsersReadLine is False:
                                        lastLine = self.recordList.readline()
                    if onQuitReplayRunFunction is not None:
                        onQuitReplayRunFunction[replayCount](parser)
                    replayCount+=1

        if onQuitRunFunction is not None:
            onQuitRunFunction(parser)

    def fileNameListGroupMatcher(self):
        """Builds the Sub Matchers for the Main Parser
        """
        return [r"(?:\||\s*)\s*(?:"
                "|".join(
                    ["%s(?:=|:)\s*(?P<%s>[0-9a-zA-Z_./\-]+)" % (fileNL.matchStr, 
                        fileNL.metaHeader + '_' + fileNL.metaNameTag + '_' + fileNL.metaName) 
                        for fileNL in self.fileDict.values()
                    ]) + ")"
                ]

    def mainFileDescription(self):
        # assemble matchers and submatchers
        return SM(name='Root',
            startReStr="",
            forwardMatch=True,
            weak=True,
            subMatchers=self.build_subMatchers()
            ) # END Root

    def build_subMatchers(self):
        return []


