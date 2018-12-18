from __future__ import print_function
from __future__ import division
from past.builtins import cmp
from builtins import zip
from builtins import range
from past.utils import old_div
import readline # optional, will allow Up/Down/History in the console
import code
import re, glob, sys, os, os.path, math, json

def goInteractive(locals):
    """Debugging function that when called stops execution and drops in an interactive loop.

Exiting the interpreter will continue execution.
call as follows:
    goInteractive(locals())
"""
    vars = globals().copy()
    vars.update(locals)
    shell = code.InteractiveConsole(vars)
    shell.interact()

def duLineGBSize(line):
    """returns the size in GB of a line outputted by du -h"""
    sizeRe=re.compile("\s*(?P<nr>[0-9.]+)(?P<unit>[kMGT]?)")
    m = sizeRe.match(line)
    if m:
        unit = {
            "": 1.0/(1024*1024*1024),
            "K": 1.0/(1024*1024),
            "M":1.0/1024,
            "G":1.0,
            "T":1024.0
        }
        return float(m.group("nr"))*unit[m.group("unit")]
    else:
        return 0.0

def totalGBSize(filePath):
    """returns the total size of the contents of a file containing du -h lines

    Makes a simple sum, which is wrong if both the directories and their contents
    are listed"""
    totS = 0.0
    with open(filePath) as f:
        while True:
            line = f.readline()
            if not line: break
            totS += duLineGBSize(line)
    return totS

def makeSplits(filePath, outSplitsF):
    groups = readDuFileToGroups(filePath)
    splits = splitGroups(groups)
    printSplits(splits)
    writeSplits(splits, outSplitsF)

def readDuFileToGroups(filePath, basePath = None):
    if basePath is None:
        basePath = os.getcwd()
    groups={}
    duLineRe = re.compile(r"\s*(?P<size>[0-9.]+[KMGT]?)\s+(?P<dir>.+)$")
    with open(filePath) as f:
        while True:
            line = f.readline()
            if not line: break
            if not line.isspace():
                m = duLineRe.match(line)
                if not m:
                    print("failed to match", repr(line))
                directory = os.path.normpath(os.path.join(basePath, m.group("dir")))
                group = os.path.dirname(directory)
                dirInGroup = groups.get(group, None)
                if dirInGroup:
                    dirInGroup[directory] = duLineGBSize(line)
                else:
                    groups[group] = {
                        directory: duLineGBSize(line)
                    }
    return groups

def splitGroups(groups, joinSmall = True):
    smallGroupsSizes = {}
    splits = []
    for group, dirInGroup in groups.items():
        groupSize = sum(dirInGroup.values())
        if joinSmall and groupSize < 0.5:
            parent = os.path.dirname(group)
            smallInParent = smallGroupsSizes.get(parent, None)
            if smallInParent:
                smallInParent[group] = groupSize
            else:
                smallGroupsSizes[parent] = {
                    group: groupSize
                }
        elif groupSize < 16:
            contents = {}
            for dirPath, size in dirInGroup.items():
                contents[os.path.basename(dirPath)] = size
            splits.append({
                "baseDir": group,
                "contents": contents,
                "size": groupSize
            })
        else:
            items = list(dirInGroup.items())
            items.sort(lambda x, y: cmp(x[1], y[1]))
            nSplits = int(math.ceil(groupSize / 14))
            splitsEls = []
            for i in range(nSplits):
                splitsEls.append({})
            while items:
                for i in range(nSplits):
                    if not items: break
                    k, v = items.pop()
                    splitsEls[i][k] = v
            for split in splitsEls:
                contents = {}
                for dirPath, size in split.items():
                    contents[os.path.basename(dirPath)] = size
                splits.append({
                        "baseDir": group,
                        "contents": contents,
                        "size": sum(split.values())
                    })
    if len(smallGroupsSizes) == 0:
        pass
    else:
        splits.extend(splitGroups(smallGroupsSizes, joinSmall = False))
    return splits

def printSplits(splits, fOut = sys.stdout):
    splits = list(splits)
    splits.sort(lambda x, y: cmp(x["size"], y["size"]))
    sizes = [x.get("size", 0) for x in splits]
    sumSize = sum(sizes)
    json.dump(splits, fOut, indent = 2)
    fOut.write("\nmin: %6.3f max: %6.3f avg: %6.3f (%7.3f for %d blocks)\n" %
               (sizes[0], sizes[-1], sumSize/len(splits), sumSize, len(splits)))

def writeSplits(splits, fOut = sys.stdout):
    for split in splits:
        fOut.write("\n")
        fOut.write(split["baseDir"])
        fOut.write("\n")
        contents = split.get("contents")
        if contents:
            for c in contents.keys():
                fOut.write(c)
                fOut.write("\n")


def filesGbSizes(globPattern):
    splitSizes = {}
    for f in glob.glob(globPattern):
        splitSizes[f] = totalGBSize(f)
    return splitSizes

def writeFilesGbSizes(globPattern):
    "prints the sizes of the files matching the glob pattern"""
    minSize = None
    maxSize = 0.0
    sumSize = 0.0
    nSizes = 0
    fileSizes = list(filesGbSizes(globPattern).items())
    fileSizes.sort(lambda x,y: cmp(x[1],y[1]))
    for f, sizeNow in fileSizes:
        print("%6.3f %s" % (sizeNow , f))
        if minSize is None or minSize > sizeNow:
            minSize = sizeNow
        if maxSize < sizeNow:
            maxSize = sizeNow
        nSizes += 1
        sumSize += sizeNow
    print()
    print("min: %6.3f max: %6.3f avg: %6.3f (%7.3f for %d blocks)" % (
        minSize, maxSize, sumSize / nSizes, sumSize, nSizes))

def writeOut(outF):
    for f in glob.glob("aflowLib3splits*"):
        outF.write("/aflowlib_data/LIB3_LIB/\n")
        for line in open(f).readlines():
            cont = transformLine(line)
            if cont:
                outF.write(cont)
                outF.write("\n")
        outF.write("\n")

lineRe=re.compile("[0-9.]+[KMGT]?\t.*/([^/]+)\n$")
lineRe=re.compile("[0-9.]+[KMGT]?\tLIB3_LIB/([^/]+)\n$")
def transformLine(line):
    m = lineRe.match(line)
    if m:
        return m.groups()[0]
    else:
        if not line.isspace():
            raise Exception("invalid line %s" % line)
        return None

def oqmdSplit(oqmdSplitFile, filesToWrite = None):
    if filesToWrite is None:
        filesToWrite = ["bryce", "perovskite", "thermoelectrics", "elements", "voorhees", "garnet", "ilmenite", "scott", "spinel", "cod", "heusler", "james"]
    # "icsd", "enum_all"
    fs = {}
    s=re.compile("/")
    r=re.compile(r"[0-9.]+[KMGT]?\s\./ftp_upload_for_uid_290/oqmd/(?P<topLevel>[^/]+).*")
    for f in filesToWrite:
        fs[f] = open("ftp_upload_for_uid_290." + f, "w")
    depth = { "binaries": 5, "icsd": 4, "bryce":5, "perovskite":5, "thermoelectrics": 5, "elements": 5, "voorhees": 4, "enum_all": 7, "garnet":5, "ilmenite":4, "scott":4, "spinel": 5, "cod": 4, "cod": 4, "heusler": 5, "james":5 }
    while True:
        line = oqmdSplitFile.readline()
        if not line: break
        m = r.match(line)
        if m:
            topLevel=m.group("topLevel")
            outF = fs.get(topLevel, None)
            if outF and depth[topLevel] == len(s.findall(line)):
                outF.write(line)

def readSizes(duFilePath):
    "return the sizes contained in a du file"
    allSizes = []
    with open(duFilePath) as fIn:
        while True:
            line = fIn.readline()
            if not line: break
            size = duLineGBSize(line)
            allSizes.append(size)
    return allSizes

def binGBSizes(allSizes):
    """returns the sizes binned in logarithmic spaced bins"""
    bins = [(1024**-2,"1K"),(10*1024**-2,"10K"),(100*1024**-2,"100K"),(1024**-1,"1M"),(10*1024**-1,"10M"), (100*1024**-1,"100M"),(1,"1G"),(10,"10G"),(100,"100G"),(1024,"1T"),(10*1024,"10T")]
    binnedValues = [ 0 for i in range(len(bins))]
    cumulativeSizes = [0.0 for i in range(len(bins))]
    for value in allSizes:
        i = 0
        while value > bins[i][0]:
            i += 1
        binnedValues[i] += 1
        cumulativeSizes[i] += value
    return { "binnedValues": binnedValues, "cumulativeValue": cumulativeSizes, "bins": bins }

def binPairValues(values, bins):
    binnedValues = [ 0 for i in range(len(bins))]
    for value, nVal in values:
        i = 0
        while value > bins[i]:
            i += 1
        binnedValues[i] += nVal
    return { "binnedValues": binnedValues, "bins": bins }

kindDist=[u'nElementalMetalsR', u'nBinaryAlloyR' ,  u'nSemiconductorsR', u'nMetalsR' , u'nRuns']

def printOut(dict):
    for k, v in dict.items():
        sys.stdout.write("%s\t%s\n" % (k, v))

bins = [3.333,10,33.33,100,333.3,1000,3333,10000, 33333, 100000]

def toBinPairs(hist):
    return [(int(x[0]), x[1]) for x in list(hist.items())]

def printHist(histDict, outF = sys.stdout):
    res=binPairValues(toBinPairs(histDict), bins)
    for b,v in zip(res["bins"], res["binnedValues"]):
        outF.write("%s\t%s\n" % (b, v))

def printTopValues(topValues, fOut = sys.stdout):
    for n,l in topValues:
        fOut.write("%s\t%s\n" %(" ".join(l), n))
