import setup_paths
import numpy as np
from contextlib import contextmanager
import logging
import json
import os
import re
import ast
from collections import namedtuple

COMMON_META_INFO_PATH = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 
    "../../../../../nomad-meta-info/meta_info/nomad_meta_info/common.nomadmetainfo.json"))

PUBLIC_META_INFO_PATH = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 
    "../../../../../nomad-meta-info/meta_info/nomad_meta_info/public.nomadmetainfo.json"))

NOTEXCEPT = re.compile(r'[a-cf-zA-CF-Z!\?,{}\[\]]')
    
def is_number(val):
        try:
            float(val)
            return True
        except ValueError:
            return False

def strcleaner(val):
        unwantedkeys = [ 
                "system", "eval", "return",
                "ctypes", "setup", "import", 
                "git", "swig", "cython"]
        for keyword in unwantedkeys:
            val = val.replace(keyword, '')
        return val

def strisinstance(val, typ):
        typlist=[]
        if isinstance(typ, str):
            typlist.append(typ)
        elif isinstance(typ, (tuple,list)):
            for ty in typ:
                typlist.append(ty)
        anytype=None
        for t in typlist:
            try:
                if("list" in t and 
                   "[" in val and 
                   "]" in val and 
                   "," in val):
                    if isinstance(literal_eval(strcleaner(val)), list):
                        anytype="list"
                        break
                    elif isinstance(literal_eval(strcleaner(val)), np.ndarray):
                        anytype="np.ndarray"
                        break
                elif("tuple" in t and 
                   "(" in val and 
                   ")" in val and 
                   "," in val):
                    if isinstance(literal_eval(strcleaner(val)), tuple):
                        anytype="tuple"
                        break
            except (TypeError,ValueError,AttributeError):
                pass
        return anytype

def literal_eval(val):
        try:
            return ast.literal_eval(val)
        except (TypeError,ValueError):
            return val

class Container(object):
    """The container class for nested data storage
    """
    def __init__(self, *args):
        self.Name = []
        for arg in args:
            if isinstance(arg, str):
                self.Name.append(arg)
        self.Active = False
        self.OpenBackend = False
        self.CloseBackend = False
        self.gid = None
        self.Containers = []
        self.References = []
        self.ReferencedFrom = []
        self.Storage = None
        self.Indent = '   '
        self.Color = None
        self.PrintOnlyActive = None
        self.PrintOnlyNames = None
        self.localDict = None
        self.currentValue = None

    def add(self, *args):
        for arg in args:
            if isinstance(arg, Container):
                self.Containers.append(arg) 
            if isinstance(arg, Storage):
                if self.Storage:
                    self.Storage(arg.__dict__)
                else:
                    self.Storage = arg
            if isinstance(arg, dict):
                if self.Storage:
                    if isinstance(self.Storage, dict):
                        self.Storage.update(arg)
                    elif isinstance(self.Storage, Storage):
                        self.Storage(arg)
                else:
                    self.Storage = Storage(arg)
            if isinstance(arg, JsonMetaInfo):
                for item in arg.jsonList:
                    self.add(item) 

    def build(self, metadata, startsection, umaskdict):
        if isinstance(metadata, JsonMetaInfo):
            attrdict = metadata.attributes(startsection)
            if attrdict:
                self.add(attrdict)
            childs = metadata.siblings(startsection)
            if startsection in umaskdict.keys():
                excludes = umaskdict[startsection]
            else:
                excludes = []
            for section in childs:
                if section not in excludes:
                    newContainer = Container(section)
                    newContainer.build(metadata, section, umaskdict)
                    if newContainer.Storage is None:
                        newContainer.add()
                    self.add(newContainer)
            return True
        else:
            return False

    def populate(self, metadata, startsection, umaskdict, updatedict):
        if isinstance(metadata, JsonMetaInfo):
            self.build(metadata, startsection, umaskdict)
            self.update(updatedict)

    def updaterefs(self, *args):
        for arg in args:
            if isinstance(arg, Container):
                self.References.extend(arg.Name) 

    def update(self, *args):
        for arg in args:
            if isinstance(arg, dict):
                if arg["startSection"]:
                    if self.Name in arg["startSection"]:
                        self.accumulateValues(self, arg)
                    else:
                        if self.Containers:
                            for module in self.Containers:
                                module.update(arg)
                else:
                    self.accumulateValues(self, arg)

    def reset(self, *args):
        for arg in args:
            if isinstance(arg, dict):
                if arg["startSection"]:
                    if self.Name in arg["startSection"]:
                        self.resetValues(self, arg)
                    else:
                        if self.Containers:
                            for module in self.Containers:
                                module.reset(arg)
                else:
                    self.resetValues(self, arg)

    def updateBackend(self, backend, startsection=None, autoopenclose=False):
        if startsection:
            if self.Name == startsection:
                if autoopenclose:
                    with self.autosection(backend, self.Name):
                        self.updateBackendValues(backend)
                else:
                    self.updateBackendValues(backend)
            else:
                if self.Containers:
                    for module in self.Containers:
                        module.updateBackend(backend, startsection, autoopenclose)
        else:
            if autoopenclose:
                with self.autosection(backend, self.Name):
                    self.updateBackendValues(backend)
            else:
                self.updateBackendValues(backend)

    @contextmanager
    def autosection(self, backend, name):
        self.gid = backend.openSection(name)
        yield self.gid
        backend.closeSection(name, self.gid)

    def opensection(self, backend, name):
        self.gid = backend.openSection(name)
        yield self.gid

    def closesection(self, backend, name):
        backend.closeSection(name, self.gid)

    def fetchAttr(self, resdict):
        for item in resdict:
            if self.Storage:
                if item in self.Storage.__dict__:
                    resdict.update({item: self.Storage.__dict__[item]})
                else:
                    if self.Containers:
                        for module in self.Containers:
                            resdict.update(module.fetchAttr(resdict))
            else:
                if self.Containers:
                    for module in self.Containers:
                        resdict.update(module.fetchAttr(resdict))
        return resdict

    def fetchAttrValue(self, itemName):
        rtnValue=None
        if self.Storage:
            if itemName in self.Storage.__dict__.keys():
                rtnValue=self.Storage.__dict__[itemName]
            else:
                if self.Containers:
                    for module in self.Containers:
                        rtnValue=module.fetchAttrValue(itemName)
                        if rtnValue is not None:
                            break
        else:
            if self.Containers:
                for module in self.Containers:
                    rtnValue=module.fetchAttrValue(itemName)
                    if rtnValue is not None:
                        break
        return rtnValue

    def updateBackendValues(self, backend):
        if self.Storage:
            self.updateBackendStorage(backend)
            self.Active = False
        if self.Containers:
            for module in self.Containers:
                module.updateBackendValues(backend)

    def accumulateValues(self, *args):
        for arg in args:
            if isinstance(arg, dict):
                if self.Storage:
                    self.accumulateDict(arg["dictionary"])
                if self.Containers:
                    for module in self.Containers:
                        module.accumulateValues(arg)
                if "activeSections" in arg:
                    if self.Name in arg["activeSections"]:
                        self.Active = True
                if "muteSections" in arg:
                    if self.Name in arg["muteSections"]:
                        self.Active = False
                        
    def resetValues(self, *args):
        for arg in args:
            if isinstance(arg, dict):
                if self.Storage:
                    self.resetAllValues()
                if self.Containers:
                    for module in self.Containers:
                        module.resetValues(arg)
                self.Active = False

    def checkUpdateValue(self, item, localdict):
        """ Updating value with the rules given in depends

         Updating values follows the following order:
          1) If 'depends' is supplied (not empty or not None), 
             the tests in the depends list will be checked in order.
             If one of the tests is successful than the one of the 
             values in 'assign' or 'value' will be updated for the item.
             Here 'assign' will assign a new value in the given string and
             'value' will assign the value of the given key item.
             (Ex. : 'assign' : 'CG' will update the item with 'CG' while 
             'value' : 'NATOM' will update the item with number of atoms 
             returned by the value of NATOM key ,which is stored in lookup dict.)
          2) If 'depends' is supplied but a lookup dictionary is not than 
             only the values of attributes in the sections can be used for test.
             The rest of the tests and assignments are updated as in case (1).
          3) If 'depends' is not supplied, subfunction is used to update value.
          4) If 'depends' and subfunction are not supplied but value of 
             MetaInfoMap is supplied, the value will be assign directly from the value 
             item.
          5) If none of the above items are supplied, this function will return None 
             to not update any values for the selected item.
        """
        # Check whether depends is supplied in the item.
        updateValue = None
        storeValue = False
        self.currentValue = updateValue
        if "prefunction" in item:
            prefunc = item.prefunction
            storeValue, updateValue, item = prefunc(item)
            self.currentValue = updateValue
        if "depends" in item:
            firstdepend = item["depends"][0]
            if "lookupdict" in item:
                needFetchVal = False
                if "test" in firstdepend:
                    storeValue, updateValue, localdict = self.checkTestsDicts(item, localdict)
                elif "assign" in firstdepend:
                    updateValue = firstdepend["assign"]
                elif "value" in firstdepend:
                    itemdepval = firstdepend["value"]
                    needFetchVal = True
                elif "store" in firstdepend:
                    itemdepval = firstdepend["store"]
                    needFetchVal = True
                    storeValue = True
                if needFetchVal:
                    if itemdepval in localdict:
                        checkval = localdict[itemdepval]
                    else:
                        accessName, checkval = self.findNameInLookupDict(itemdepval, item.lookupdict)
                        localdict.update({itemdepval : checkval})
                    updateValue = checkval
            else:
                needFetchVal = False
                if "test" in firstdepend:
                    storeValue, updateValue, localdict = self.checkTestsAttr(item, localdict)
                elif "assign" in firstdepend:
                    updateValue = firstdepend["assign"]
                elif "value" in firstdepend:
                    itemdepval = firstdepend["value"]
                    needFetchVal = True
                elif "store" in firstdepend:
                    itemdepval = firstdepend["store"]
                    needFetchVal = True
                    storeValue = True
                if needFetchVal:
                    if itemdepval in localdict:
                        checkval = localdict[itemdepval]
                    else:
                        attrdict = {deptest[0] : ''}
                        attrdict = self.fetchAttr(attrdict)
                        localdict.update(attrdict)
                        checkval = attrdict[deptest[0]]
                    updateValue = checkval
        self.currentValue = updateValue
        if "subfunction" in item:
            subfunc = item.subfunction["function"]
            if "supportDict" in item.subfunction:
                supDict = item.subfunction["supportDict"]
                storeValue, updateValue, item = subfunc(supDict, item)
            else:
                storeValue, updateValue, item = subfunc(item)
            self.currentValue = updateValue
        if "value" in item:
            updateValue = item['value']
            self.currentValue = updateValue
        if("depends" in item or "value" in item or "subfunction" in item):
            pass
        elif "defaultValue" in item:
            updateValue = item['defaultValue']
            self.currentValue = updateValue
        if "postfunction" in item:
            postfunc = item.postfunction
            storeValue, updateValue, item = postfunc(item)
            self.currentValue = updateValue
        if "valtype" in item:
            if updateValue is not None:
                updateValue = self.convertToNumber(updateValue, item["valtype"])
        if("unit" in item and "unitdict" in item):
            if updateValue is not None:
                updateValue = self.convertUnits(updateValue, item["unit"], item["unitdict"])
        return storeValue, updateValue, localdict

    def convertToNumber(self, updateValue, valtype):
        acceptvals = ["float", "int", "list", "tuple", 
                      "np.asarray", "np.array", 
                      "set", "boolean", "bytes"]
        if valtype in acceptvals:
            if(isinstance(updateValue, list) or isinstance(updateValue, tuple)):
                try:
                    newUpdateValue = [eval(
                        valtype+"("+literal_eval(str(ival))+")"
                        ) for ival in updateValue]
                except (TypeError,ValueError):
                    newUpdateValue = updateValue
            elif isinstance(updateValue, np.ndarray):
                try:
                    newUpdateValue = np.asarray([eval(
                        valtype+"("+literal_eval(str(ival))+")"
                        ) for ival in updateValue])
                except (TypeError,ValueError):
                    newUpdateValue = updateValue
            elif is_number(updateValue):
                newUpdateValue = eval(valtype+"("+str(updateValue)+")")
            else:
                newUpdateValue = updateValue
        else:
            # I hope you know what you are doing
            try:
                newUpdateValue = float(updateValue)
            except (TypeError,ValueError):
                newUpdateValue = updateValue
        return newUpdateValue
    
    def convertUnits(self, updateValue, unit, unitdict):
        if(isinstance(updateValue, list) or isinstance(updateValue, tuple)):
            updateValue = [float(ival) * self.unitConverter(
                unit, unitdict
                ) for ival in updateValue]
        elif isinstance(updateValue, np.ndarray):
            updateValue = updateValue * self.unitConverter(
                unit, unitdict)
        elif is_number(updateValue):
            updateValue = self.convertToNumber(updateValue, "float")
            if updateValue:
                updateValue = updateValue * self.unitConverter(
                        unit, unitdict)
        elif isinstance(updateValue, str):
            # I hope you know what you are doing
            try:
                newUpdateVal = strcleaner(updateValue)
                newUpdateVal = NOTEXCEPT.sub('', newUpdateVal)
                updateValue = float(newUpdateVal) * self.unitConverter(
                        unit, unitdict)
            except (TypeError,ValueError):
                pass
        return updateValue

    def unitConverter(self, unit, unitdict):
        """ Unit converter using definitions of units explicitly

            The unit names are converted to numbers and the resulting
            expression will be evaluated by python.
            Ex.: unit = 'electron-volt/Angstrom^3' 
                 will be converted to
                 unit = '1.602176565e-19*1.0/1.0e-10**3'
                 factor = eval(unit) = 160.2176565e+9 Joule/meter^3 (Pascal)
                 160.2176565e+9 Pascal = 160.2176565 GPa
                 in SI units and the result will be calculated as follows:
                 output_value = input_value * factor
        """
        newunit = unit.lower()
        newunit = newunit.replace('-','*').replace(' ', '*').replace('^', "**")
        for key,value in unitdict.items():
            newunit = newunit.replace(str(key), str(value))
        newunit = NOTEXCEPT.sub('', newunit)
        try:
            return float(eval(newunit))
        except (ValueError,TypeError):
            return None

    def checkTestsDicts(self, item, localdict):
        for depdict in item["depends"]:
            deptests = depdict["test"]
            depmeet = 0
            for deptest in deptests:
                if deptest[0] in localdict:
                    checkval = localdict[deptest[0]]
                else:
                    accessName, checkval = self.findNameInLookupDict(deptest[0], item.lookupdict)
                    localdict.update({deptest[0] : checkval})
                if(('<' in deptest[1] or   # In Python 3, different type comparisons
                    '>' in deptest[1]) and # are removed. Therefore, < and > comparisons
                    (checkval is None)):   # with a None value generates TypeError
                    pass
                else:
                    if isinstance(checkval, str):
                        if eval('"' + str(checkval) + '"' + deptest[1]):
                            depmeet += 1
                    else:
                        if eval(str(checkval) + deptest[1]):
                            depmeet += 1
                if depmeet == len(deptests):
                    storeValue = False
                    if 'assign' in depdict:
                        return storeValue, depdict['assign'], localdict
                    elif 'value' in depdict:
                        if depdict['value'] in localdict:
                            checkval = localdict[depdict['value']]
                        else:
                            accessName, checkval = self.findNameInLookupDict(depdict['value'], 
                                    item.lookupdict)
                            localdict.update({depdict['value'] : checkval})
                        return storeValue, checkval, localdict
                    elif "store" in depdict:
                        itemdepval = depdict["store"]
                        storeValue = True
                        if itemdepval in localdict:
                            checkval = localdict[itemdepval]
                        else:
                            accessName, checkval = self.findNameInLookupDict(itemdepval, item.lookupdict)
                            localdict.update({itemdepval : checkval})
                        print("PRINTING: metainfo, store:",storeValue, checkval)
                        return storeValue, checkval, localdict
        return False, None, localdict

    def checkTestsAttr(self, item, localdict):
        for depdict in item["depends"]:
            #depdict = item["depends"][tests]
            for deptests in depdict["test"]:
                depmeet = 0
                for deptest in deptests:
                    if deptest[0] in localdict:
                        checkval = localdict[deptest[0]]
                    else:
                        attrdict = {deptest[0] : ''}
                        attrdict = self.fetchAttr(attrdict)
                        localdict.update(attrdict)
                        checkval = attrdict[deptest[0]]
                    if eval(str(checkval) + deptest[1]):
                        depmeet += 1
                if depmeet == len(deptests):
                    if 'assign' in depdict:
                        return depdict['assign'], localdict
                    elif 'value' in depdict:
                        if depdict['value'] in localdict:
                            checkval = localdict[depdict['value']]
                        else:
                            attrdict = {depdict['value'] : ''}
                            attrdict = self.fetchAttr(attrdict)
                            localdict.update(attrdict)
                            checkval = attrdict[deptest[0]]
                        return checkval, localdict
        return None, localdict

    def findNameInLookupDict(self, metaname, lookupdict):
        for item in lookupdict:
            itemMap = lookupdict[item]
            #if metaname in itemMap.metaName:
            if metaname in itemMap.matchStr:
                return item, itemMap.value
        return None, None

    def updateBackendStorage(self, backend):
        for itemk in self.Storage.__dict__:
            if(self.Storage.__dict__[itemk]["act"] or 
                self.Storage.__dict__[itemk]["stor"]):
                self.Storage.__dict__[itemk]["act"] = False
                self.Storage.__dict__[itemk]["stor"] = False
                value = self.Storage.__dict__[itemk]["val"]
                if isinstance(value, np.ndarray):
                    backend.addArrayValues(itemk, value)
                elif isinstance(value, (list, tuple)):
                    backend.addArrayValues(itemk, np.asarray(value))
                elif value is None:
                    pass
                else:
                    backend.addValue(itemk, value)

    def accumulateDict(self, checkDict):
        localdict = {}
        for itemk in checkDict:
            itemv = checkDict[itemk]
            storeValue, updateValue, localdict = self.checkUpdateValue(itemv, localdict)
            if updateValue is not None:
                itemkinlist = True if itemk in self.Storage.__dict__ else False
                if itemkinlist is False:
                    if (itemk.startswith('x_') and 
                        itemv.activeSections): 
                        attr_update=False
                        for actSect in itemv.activeSections:
                            if actSect in self.Name:
                                attr_update=True
                        if attr_update:
                            attr_act = False
                            attrvalues = {
                                'act' : attr_act, 
                                'val' : None, 
                                'stor': False, 
                                'kind': None, 
                                'dtyp': "C",
                                'unit': None,
                                'size': [],
                                'refs': None
                                }
                            self.Storage.__dict__.update({itemk : attrvalues})
                if itemk in self.Storage.__dict__:
                    if storeValue:
                        #If we need to store the updated values
                        if self.Storage.__dict__[itemk]["val"] is None:
                            #If not initialized, initialize with a list to store
                            if isinstance(updateValue, str) and float(updateValue):
                                self.Storage.__dict__[itemk]["val"] = [float(updateValue)]
                                self.Storage.__dict__[itemk]["stor"] = True
                            else:
                                self.Storage.__dict__[itemk]["val"] = [updateValue]
                                self.Storage.__dict__[itemk]["stor"] = True
                        else:
                            #Append to the stored list if there is a list
                            if type(self.Storage.__dict__[itemk]["val"]) is list:
                                #self.Storage.__dict__[itemk]["val"].append(updateValue)
                                if isinstance(updateValue, str) and float(updateValue):
                                    self.Storage.__dict__[itemk]["val"].append(float(updateValue))
                                    self.Storage.__dict__[itemk]["stor"] = True
                                else:
                                    self.Storage.__dict__[itemk]["val"].append(updateValue)
                                    self.Storage.__dict__[itemk]["stor"] = True
                            else:
                                #Convert the prevoius update to list and append update
                                preValue = self.Storage.__dict__[itemk]["val"]
                                self.Storage.__dict__[itemk]["val"] = [preValue, updateValue]
                                self.Storage.__dict__[itemk]["stor"] = True
                    else:
                        #No need to store, assign the updated value
                        self.Storage.__dict__[itemk]["val"] = updateValue
                    if itemv.activeInfo:
                        if (itemk.startswith('x_') and 
                            itemv.activeSections): 
                            attr_update=False
                            for actSect in itemv.activeSections:
                                if actSect in self.Name:
                                    attr_update=True
                            if attr_update:
                                self.Storage.__dict__[itemk]["act"] = True
                        else:
                            self.Storage.__dict__[itemk]["act"] = True
                    else:
                        self.Storage.__dict__[itemk]["act"] = False
                self.Active = True
                if "valueSize" in itemv:
                    if "sizeMetaName" in itemv:
                        self.Storage.__dict__[itemv["sizeMetaName"]] = itemv["valueSize"]
                if "unitconverter" in itemv:
                    newValue = itemv["unitconverter"](self, itemv)
                    self.Storage.__dict__[itemk["val"]] = newvalue

    def resetAllValues(self):
        for itemk in self.Storage.__dict__:
            self.Storage.__dict__[itemk]["val"] = None
            self.Storage.__dict__[itemk]["stor"] = False
            self.Storage.__dict__[itemk]["act"] = False
            self.Active = False

    def __str__(self, caller=None, decorate='', color=None, printactive=None, onlynames=None):
        string = ''
        if onlynames is None:
            if self.PrintOnlyNames:
                onlynames = self.PrintOnlyNames
        if printactive is None:
            if self.PrintOnlyActive:
                printactive = self.PrintOnlyActive
        if color:
            color = int(color) + 1
        else:
            if self.Color:
                color = int(self.Color)
        printok = False
        if printactive:
            if self.Active:
                printok = True
        else:
            printok = True
        if caller:
            if printok:
                if color:
                    string = '%s\033[9%sm`-->[' % (decorate + self.Indent, str(color%6 + 1))
                    string = string + ','.join(['%s' % (name) for name in self.Name]) + ']\n\033[0m'
                else:
                    string = '%s`-->[' % (decorate + self.Indent)
                    string = string + ','.join(['%s' % (name) for name in self.Name]) + ']\n'
        else:
            if printok:
                if color:
                    string = '%s\033[9%sm-->[' % (decorate + self.Indent, str(color%6 + 1)) 
                    string = string + ','.join(['%s' % (name) for name in self.Name]) + ']\n\033[0m'
                else:
                    string = '%s-->[' % (decorate + self.Indent) 
                    string = string + ','.join(['%s' % (name) for name in self.Name]) + ']\n'
        if printok:
            if color:
                string = string + '%s\033[9%sm|\033[0m   `.\n' % (decorate + self.Indent, str(color%6 + 1)) 
            else:
                string = string + '%s|   `.\n' % (decorate + self.Indent) 
        if self.Storage:
            for key in self.Storage.__dict__:
                printattrok = False
                if printactive:
                    if self.Storage.__dict__[key]["act"]:
                        printattrok = True
                else:
                    printattrok = True
                if printattrok and printok:
                    if color:
                        if onlynames:
                            string = string + '%s\033[9%sm|\033[0m    |__.%s\n' % (decorate 
                                    + self.Indent, str(color%6 + 1), key)
                        else:
                            string = string + '%s\033[9%sm|\033[0m    |__.%s : Active=%s Value=%s\n' % (decorate 
                                    + self.Indent, str(color%6 + 1), key, self.Storage.__dict__[key]["act"], 
                                    self.Storage.__dict__[key]["val"])
                    else:
                        if onlynames:
                            string = string + '%s|    |__.%s\n' % (decorate + 
                                    self.Indent, key)
                        else:
                            string = string + '%s|    |__.%s : Active=%s Value=%s\n' % (decorate + 
                                    self.Indent, key, self.Storage.__dict__[key]["act"], 
                                    self.Storage.__dict__[key]["val"])
            if color:
                string = string + '%s\033[9%sm|\033[0m\n' % (decorate + self.Indent, str(color%6 + 1))
            else:
                string = string + '%s|\n' % (decorate + self.Indent)
        if self.Containers:
            for module in self.Containers:
                if color:
                    string = string + '%s\033[9%sm|\033[0m\n' % (decorate + self.Indent, 
                            str(color%6 + 1)) + module.__str__(self.Name, 
                            '%s\033[9%sm|\033[0m' % (decorate + self.Indent, str(color%6 + 1)), 
                            color, printactive, onlynames)
                else:
                    string = string + '%s|\n' % (decorate + self.Indent) + module.__str__(self.Name, 
                            '%s|' % (decorate + self.Indent), 
                            color, printactive, onlynames)
        return string

class Storage(dict):
    """ Storage for meta info document types.
        Sections are build by Container class
    """
    def __init__(self, *args, **kwargs):
        super(Storage, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    self[k] = v

        if kwargs:
            for k, v in kwargs.items():
                self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(Storage, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(Storage, self).__delitem__(key)
        del self.__dict__[key]


class JsonMetaInfo(object):
    """ Json file loader for meta info data of NOMAD.
        Loads data and extracts values of items 
        with specified superNames
    """
    def __init__(self, *args):
        self.jsonList = None
        for filepath in args:
            try:
                with open(filepath, encoding="utf-8") as f:
                    jsonDict = json.load(f)
            except:
                logging.exception("Error while loading file %s" % filepath)
                raise
            typeStr = jsonDict.get("type","nomad_meta_info_1_0")
            typeRe = re.compile(r"nomad_meta_info_(?P<major>[0-9]+)_(?P<minor>[0-9]+)$")
            m = typeRe.match(typeStr)
            if not m:
                raise Exception("unexpected type '%s', expected nomad_meta_info_1_0" % typeStr)
            newJsonList = jsonDict.get("metaInfos",[])
            if self.jsonList:
                self.jsonList = self.jsonList + newJsonList
            else:
                self.jsonList = newJsonList

    def attributes(self, sectionname):
        attributes = {}
        for item in self.jsonList:
            superlist = item['superNames']
            itemname = item['name']
            try:
                kindname = item['kindStr']
            except:
                kindname = []
            try:
                dtyp = item['dtypeStr']
            except:
                dtyp = []
            try:
                size = item['shape']
            except:
                size = []
            try:
                unit = item['units']
            except:
                unit = []
            try:
                refs = item['referencedSections']
            except:
                refs = []
            if ('type_section' in kindname or 
                sectionname not in superlist):
                continue
            attrvalues = {
                    'act' : False, 
                    'stor': False, 
                    'val' : None, 
                    'kind': kindname, 
                    'dtyp': dtyp,
                    'unit': unit,
                    'size': size,
                    'refs': refs
                    }
            attributes.update({itemname: attrvalues})
        return attributes

    def siblings(self, sectionname):
        siblings = []
        searchList = []
        nameList = []
        containsList = []
        for item in self.jsonList:
            superlist = item['superNames']
            itemname = item['name']
            try:
                kindname = item['kindStr']
            except:
                kindname = []
            if ('type_section' in kindname or 
                'type_abstract_document_content' in kindname): 
                if sectionname in superlist:
                    searchList.append(itemname)
                if itemname not in nameList:
                    nameList.append(itemname)
        for name in searchList:
            if (set([name]) not in set(nameList) and 
                self.isparent(name)):
                siblings.append(name)
        return siblings

    def isparent(self, itemname):
        haschild = False
        for item in self.jsonList:
            if itemname in item['superNames']:
                haschild = True
        return haschild

    def rootsections(self):
        rootname = []
        searchList = []
        nameList = []
        for item in self.jsonList:
            superlist = item['superNames']
            itemname = item['name']
            try:
                kindname = item['kindStr']
            except:
                kindname = []
            if 'type_section' not in kindname:
                continue
            if not superlist:
                searchList.append(itemname)
            if itemname not in nameList:
               nameList.append(itemname)
        for name in searchList:
            if ('section' in name and set([name]) not in set(nameList)):
                rootname.append(name)
        return rootname

    def fetchdict(self, itemname, pattern):
        resDict = {}
        for item in self.jsonList:
            val = dict(item)
            itemProperty = item[itemname]
            if pattern in itemProperty:
                resDict.update({item["name"]: val})
        return resDict

if __name__ == "__main__":
    run = Container('section_run')
    exclude_dict = { 
            'section_run' : [
            'section_processor_info', 
            'section_processor_log', 
            'section_springer_material',
            'section_repository_info'
            ]}

    jsonmetadata = JsonMetaInfo(COMMON_META_INFO_PATH, PUBLIC_META_INFO_PATH)

    updateDict = {
            'startSection' : [['section_topology']],
            'muteSections' : [['section_interaction']],
            'dictionary' : {
                'molecule_constraint_atoms' : {'depends' : [[]], 'assign' : 100},
                'interaction_atoms' : {'depends' : [[]], 'assign' : 10},
                'topology_force_field_name' : {'depends' : [[]], 'assign' : "ReaxFF"}
                }
            }
    
    run.populate(jsonmetadata, 'section_run', exclude_dict, updateDict)
    run.Color = 4
    for container in run.Containers:
        if 'section_topology' in container.Name:
            select = container
    select.Color = 4
    #select.PrintOnlyActive = 1
    run.PrintOnlyNames = 1
    print(run)


