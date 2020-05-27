import setup_paths

def get_unitDict(keyname):
    """ Unit dictionary for convertions

        Unit names will be converted to values.
        When defining units in translator dictionary, 
        the unit names in the dictionary should be used.
        The unit convertion values are written for SI units.
        If you would like to change it, just add another key 
        to the dictionary and change the key at parser base.
        Usage:
           Natively support python language in definitions. 
           You can use any python math operator and function. 
           Moreover, you can use space or - for multiplication 
           and ^ for power of values. 
        Example:
            kilogram/meter^2 can be written as
            kilo-gram/meter^2 or kilo-gram/meter**2
            and will be calculated as
            kilo*gram/meter**2

    """
    unitDict = {
        "si" : {
            "meter"          : "1.0",
            "kilo"           : "1.0e3",
            "gram"           : "1.0e-3",
            "second"         : "1.0",
            "joule"          : "1.0",
            "newton"         : "1.0",
            "kelvin"         : "1.0",
            "pascal"         : "1.0",
            "coulomb"        : "1.0",
            "volt"           : "1.0",
            "centi"          : "1.0e-2",
            "milli"          : "1.0e-3",
            "micro"          : "1.0e-6",
            "nano"           : "1.0e-9",
            "pico"           : "1.0e-12",
            "femto"          : "1.0e-15",
            "atto"           : "1.0e-18",
            "erg"            : "1.0e-7",
            "dyne"           : "1.0e-5",
            "bar"            : "1.0e-1",
            "angstrom"       : "1.0e-10",
            "kcal"           : "4184.096739614824",
            "mol"            : "0.602213737699784e24",
            "atmosphere"     : "1.01325e5",
            "electron"       : "1.602176565e-19",
            "atomicmassunit" : "1.66054e-27",
            "amu"            : "1.66054e-27",
            "bohr"           : "5.29177249e-11",
            "hartree"        : "4.35974e-18",
            "pascal"         : "1.0",
            "akmatime"       : "0.048888e-12",
            },
        "amber" : {
            "time"           : "pico*second",
            "length"         : "nano*meter",
            "temperature"    : "Kelvin",
            "energy"         : "kcal/mol",
            },
        "gromacs" : {
            "time"           : "pico*second",
            "length"         : "nano*meter",
            "temperature"    : "Kelvin",
            "energy"         : "kcal/mol",
            },
        "mdtraj" : {
            "time"           : "pico*second",
            "length"         : "nano*meter",
            "angle"          : "degree",
            },
        "mdanalysis" : {
            "time"           : "pico*second",
            "length"         : "angstrom",
            "energy"         : "kilo*joule/mol",
            "charge"         : "electron",
            "force"          : "kilo*joule/(mol*angstrom)",
            "velocity"       : "angstrom/(pico*second)",
            },
        "ase" : {
            "time"           : "pico*second",
            "length"         : "angstrom",
            "energy"         : "kilo*joule/mol",
            "charge"         : "electron",
            "force"          : "kilo*joule/(mol*angstrom)",
            "velocity"       : "angstrom/(pico*second)",
            },
        #AKMA units: http://www.esi.umontreal.ca/accelrys/life/insight2000.1/charmm_principles/Ch01_intro.FM5.html
        "akma" : {
            "time"           : "0.048888e-12*second",
            "outputtime"     : "pico*second",
            "length"         : "angstrom",
            "energy"         : "kcal/mol",
            "charge"         : "electron",
            "mass"           : "1.661e-27*kilo-gram",
            "force"          : "kcal/(mol*angstrom)",
            "velocity"       : "angstrom/akmatime",
            },
        }
    if keyname:
        resDict = unitDict[keyname]
    else:
        resDict = unitDict["si"]
    return resDict

def metaNameConverter(keyName):
    newName = keyName.lower().replace(" ", "").replace("-", "")
    newName = newName.replace("(", "").replace(")", "")
    newName = newName.replace("[", "").replace("]", "")
    newName = newName.replace(",", "").replace(".", "")
    newName = newName.replace("\\", "").replace("/", "")
    newName = newName.replace("'", "").replace(":", "")
    return newName

def metaNameConverter_UnderscoreSpaceDash(keyName):
    newName = keyName.lower()
    newName = ' '.join(newName.split())
    newName = newName.replace(" ", "_").replace("-", "_")
    newName = metaNameConverter(newName)
    return newName

class MetaInfoMap(dict):
    """Map cache values to meta info
    """
    activeInfo=False # The key will be activated when it is matched
    infoPurpose=None # Same as in FileMapDict
    defaultValue=None # Default value of this key
    changeTags=None # Also set the keys in the dict to given values with the match
    nameTranslate=None # translate the key to meta name using this function
    matchStr=None # this string will be matched (can be altered after initialized)
    matchNames=None # Control list/tuple/dict/str for matchStr (See getDict_MetaStrInDict)
    nextMatch=None # Do not alter the match of item until the nextMatch string is matched
    concatMatch=None # If there is an other match in list (Ex.:nextMatch) concatanate it 
    replaceTag=None # Replace the original tag with this value but still use the original matchStr
    matchWith=None # match control: Next word (NW), Previous word (PW), End of line(EOL), Until delimeter (UD)
    removeText=None # remove the given text to generate meta info name
    alsoMatch=None # also match the strings in this list with the given key in same line (any where)
    metaHeader=None # the general meta name header (See xxxDictionary for usage)
    metaName=None # this will be the meta info name that will be passed to backend and written output
    metaNameTag=None # the meta info name body 
    metaInfoType=None # Float, string or other that have to be defined in meta info core
    value=None # The final value that will be passed to backend and written to output
    valueSize=None # The size of the value that will be pased to backend and written to output
    sizeMetaName=None 
    depends=None # list/dict for altering values in this or other dictionaries
    lookupdict=None # the lookup dict for values in depends (can also be itself)
    subfunction=None # an extra update function before passing value to backend or write output
    activeSections=None # The sections that this meta info key will be active to pass backend
    autoSections=False

    def __init__(self, *args, **kwargs):
        super(MetaInfoMap, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    if k in self:
                        self[k] = v 
        if kwargs:
            for k, v in kwargs.items():
                if k in self:
                    self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(MetaInfoMap, self).__setitem__(key, value)
        self.__dict__.update({key: value})

class FileInfoMap(dict):
    """Map cache values to meta info
    """
    activeInfo=False
    infoPurpose=None
    fileName=None
    fileFormat=None 
    fileSupplied=False
    fileHolder=None
    nameTranslate=None 
    matchStr=None 
    matchNames=None 
    metaHeader=None 
    metaName=None 
    metaNameTag=None
    metaInfoType=None
    value=None 
    valueSize=None 
    sizeMetaName=None 
    depends=None
    lookupdict=None
    subfunction=None
    activeSections=None

    def __init__(self, *args, **kwargs):
        super(FileInfoMap, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    if k in self:
                        self[k] = v 
        if kwargs:
            for k, v in kwargs.items():
                if k in self:
                    self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(FileInfoMap, self).__setitem__(key, value)
        self.__dict__.update({key: value})

class MapDictionary(dict):
    """
    Modified from the reference source below:
    https://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary
    Example:
    m = MapDictionary({'Name': 'mdtraj'}, format='.mdcrd', found=True, list=['Value'])
    """
    def __init__(self, *args, **kwargs):
        super(MapDictionary, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    if (isinstance(v, FileInfoMap) or
                        isinstance(v, MetaInfoMap)):
                        if v.replaceTag:
                            k_new = v.replaceTag
                        else:
                            k_new = k
                        if v.nameTranslate:
                            v.metaName = v.nameTranslate(k_new)
                        else:
                            v.metaName = k_new
                    v.matchStr = k
                    metaStr = ''
                    if v.metaHeader:
                        metaStr = metaStr + v.metaHeader + '_'
                    if v.metaNameTag:
                        metaStr = metaStr + v.metaNameTag + '_'
                    metaStr = metaStr + v.metaName
                    self[metaStr] = v
                    if metaStr != k:
                        self.pop(k, None)

        if kwargs:
            for k, v in kwargs.items():
                if (isinstance(v, FileInfoMap) or
                    isinstance(v, MetaInfoMap)):
                    if v.metaTranslate:
                        v.metaName = v.nameTranslate(k)
                    else:
                        v.metaName = k
                v.matchStr = k
                metaStr = ''
                if v.metaHeader:
                    metaStr = metaStr + v.metaHeader + '_'
                if v.metaNameTag:
                    metaStr = metaStr + v.metaNameTag + '_'
                metaStr = metaStr + v.metaName
                self[metaStr] = v
                if metaStr != k:
                    self.pop(k, None)

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(MapDictionary, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(MapDictionary, self).__delitem__(key)
        del self.__dict__[key]

    def get_keys(self):
        return [val.metaName for val in self.__dict__.values()]

def getList_MetaStrInDict(sourceDict):
    """Returns a list that includes all meta name 
       strings for the given dictionary.
       Meta name strings are not actual meta names but 
       used as the keywords in the parsing.
    """
    return [sourceDict[item].matchStr for item in sourceDict]

def copyMetaDictToDict(sourceObj, sourceDict, sourceType=None, targetObj=None, 
                       targetDict=None, targetType=None, copyDict=None):
    """ Copy values of a source dictionary to a target one 
        by matching the keys in bith dictionaries.
        Dictionary types can be python or SmartParser

        Inputs:
            sourceObj  (object) : python object that holds sourceDict
            sourceDict (string) : the name of the source dictionary 
            sourceType (string) : type can be 'standard' or 'smartparser'
            targetObj  (object) : python object that holds targetDict
                                  if None, sourceObj will be used
            targetDict (string) : the name of the target dictionary 
                                  if None, sourceDict will be used
            targetType (string) : type can be 'standard' or 'smartparser'
                                  if targetDict is None, this equals to sourceType
            copyDict (dict)     : includes items for mapping 
                                  keys from source to target
                                  Ex. : {"ener0" : "Energy", "step" : "Timestep"}

            If target inputs are None, the copy will work on the same dictionary.

            Outputs: [Success, Status]
            Success :  True if any copy else False
            Status  : -2 if both dictionaries are not accessible
                      -1 if only target dict is not accessible
                       0 can access dict but without any copy
                       1 at least one copy is made from sourceDict to targetDict
                       N number of N copies were made from sourceDict to targetDict
    """
    rtn = [False,-2]
    skip = False
    sDict = None
    tDict = None
    if targetDict is None and copyList is None:
        skip = True
    if skip is False:
        if sourceType is None:
            sourceType = 'smartparser'
        if targetObj is None:
            targetObj = sourceObj
        if targetDict is None:
            targetDict = sourceDict
            targetType = sourceType
        if targetType is None:
            targetType = sourceType
        sDict = getattr(sourceObj, sourceDict)
        tDict = getattr(targetObj, targetDict)
        if sDict is not None:
            rtn[1] += 1
        if tDict is not None:
            rtn[1] += 1
        sp = False
        if 'smartparser' in sourceType:
            sp = True
        tp = False
        if 'smartparser' in targetType:
            tp = True
        for skey, sval in sDict.items():
            sk = sval.matchStr if sp else skey
            if copyDict:
                if sk in copyDict:
                    target = copyDict[sk]
                    for tkey, tval in tDict.items():
                        tk = tval.matchStr if tp else tkey
                        if target == tk:
                            rtn[0] = True
                            rtn[1] += 1
                            if tp:
                                tDict[tkey] = sval.value if sp else sval
                            else:
                                tDict[tkey] = sval.value if sp else sval
            else:
                for tkey, tval in tDict.items():
                    tk = tval.matchStr if tp else tkey
                    if sk == tk:
                        rtn[0] = True
                        rtn[1] += 1
                        if tp:
                            tDict[tkey] = sval.value if sp else sval
                        else:
                            tDict[tkey] = sval.value if sp else sval
                        break
        if rtn[1]>0:
            setattr(targetObj, tDict, targetDict)
    return rtn

def setMetaStrInDict(objectName, sourceDict, mStr, mVal, matchExact=False, 
        caseSensitive=False, multipleMatch=True, matchFirst=None, matchLast=None):
    """Find mStr in a given sourceDict and set it to mVal.
       If the value is set in a dictionary item, activeInfo of the 
           item is also set to True
     
       Inputs: 
           objectName (String) : Python object that holds this attribute 
                                 Ex. : self
           sourceDict (String) : The name string of the SmartParser dictionary
                                 Ex. : 'myDict'  (this will be self.myDict)
           mStr (String) : The string to be matched with the matchStr in dict.
           mVal (Any type that is accepted by metaInfoStorage) : The value
       Output 
           A list with [Success, MatchKey, Status]:
               Success:  True : there is at least one match
                         False: there is an error
               MatchKey       : a list with the matched 
                                  keys in the dict
               Status :  -1   : there is no dictionary 
                                  with the given name
                          0   : there is dict but there is no 
                                  item that matchs with mStr
                          1   : there is a match in the dict 
                                  with a successful value set
                          2-N : there are multiple matches 
                                  with successful value settings
                                  (Only available if multipleMatch is True)
    """
    rtnList = [False, None, -1]
    theDict = getattr(objectName, sourceDict)
    if theDict is not None:
        rtnList[2] = 0
        matchKeys = []
        if matchFirst is not None:
            if int(matchFirst)<1:
                matchFirst = None
            else:
                macthExact = False
        if matchLast is not None:
            if int(matchLast)<1:
                matchLast = None
            else:
                macthExact = False
        for k,v in theDict.items(): 
            assign = False
            if matchExact is True:
                if caseSensitive:
                    if mStr == v.matchStr:
                        assign = True
                else:
                    if mStr.upper() == v.matchStr.upper():
                        assign = True
            elif matchFirst is not None:
                mFirst = int(matchFirst)
                if mFirst >= len(v.matchStr):
                    mFirst=len(v.matchStr)
                vStr = v.matchStr[0:mFirst]
                if caseSensitive:
                    if mStr in vStr:
                        assign = True
                else:
                    if mStr.upper() in vStr.upper():
                        assign = True
            elif matchLast is not None:
                mLast = int(matchLast)
                if mLast >= len(v.matchStr):
                    mLast=len(v.matchStr)
                vStr = v.matchStr[-mLast:-1]
                if caseSensitive:
                    if mStr in vStr:
                        assign = True
                else:
                    if mStr.upper() in vStr.upper():
                        assign = True
            else:
                if caseSensitive:
                    if mStr in v.matchStr:
                        assign = True
                else:
                    if mStr.upper() in v.matchStr.upper():
                        assign = True
            if assign:
                theDict[k].value = mVal
                theDict[k].activeInfo = True
                matchKeys.append(k)
                if multipleMatch is False:
                    break
        if matchKeys:
            rtnList[0] = True
            rtnList[1] = matchKeys
            rtnList[2] = len(matchKeys)
            setattr(objectName, sourceDict, theDict)
    return rtnList

def isMetaStrInDict(nameStr, sourceDict):
    """Returns a list that includes all meta name 
       strings for the given dictionary.
       Meta name strings are not actual meta names but 
       used as the keywords in the parsing.
    """
    val = None
    for k,v in sourceDict.items(): 
        if nameStr in v.matchStr:
            val = k
            break
    return val
    #return [k for k,v in sourceDict.items() if nameStr in v.matchStr][0]

def getDict_MetaStrInDict(sourceDict, nameList=None):
    """Returns a dict that includes all meta name 
       strings and corresponding values for the given dictionary.
       Meta name strings are not actual meta names but 
       used as the keywords in the parsing.
    """
    newDict = {}
    if nameList is not None:
        if isinstance(nameList, (list, tuple)):
            for key in nameList:
                if key in sourceDict:
                    newDict.update({sourceDict[key].matchStr : key}) 
        elif isinstance(nameList, dict):
            for key, value in nameList:
                if value in sourceDict:
                    newDict.update({key : value}) 
        elif isinstance(nameList, str):
            for key, value in sourceDict.items():
                if nameList in sourceDict[key]:
                    if isinstance(sourceDict[key][nameList], (list, tuple)):
                        for item in sourceDict[key][nameList]:
                            newDict.update({item : key}) 
                    elif isinstance(sourceDict[key][nameList], str):
                        newDict.update({sourceDict[key][nameList] : key}) 
                    else:
                        newDict.update({sourceDict[key].matchStr : key}) 
                else:
                    newDict.update({sourceDict[key].matchStr : key}) 
    else:
        for key, value in sourceDict.items():
            newDict.update({sourceDict[key].matchStr : key}) 
    return newDict


