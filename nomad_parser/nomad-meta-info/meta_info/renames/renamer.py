#!/bin/env python

import os
import tempfile
import re
import logging
from io import open
from io import StringIO

def loadRenamesFile(renamesPath):
    """Loads the content of a rename file"""
    renamesRe = re.compile(
        r"\s*(?P<old>[a-zA-Z0-9_]+)\s*->\s*(?P<new>[a-zA-Z0-9_]+)\s*$")
    renames = {}
    with open(renamesPath, encoding='utf-8') as f:
        while True:
            line = f.readline()
            if not line:
                break
            m = renamesRe.match(line)
            if m:
                renames[m.group("old")] = m.group("new")
            elif not line.isspace():
                logging.warn("Unexpected line %r in %s", line, renamesPath)
    return renames

def renamesSearchRe(renames):
    """creates a regular expression that matches the words that should be renamed"""
    res = StringIO()
    res.write(r"\b(")
    first = True
    for k in renames.keys():
        if not first:
            res.write("|")
        else:
            first = False
        res.write(k)
    res.write(r")\b")
    renameReStr = res.getvalue()
    return re.compile(renameReStr)

def replaceInFile(filePath, replacements):
    """performs the replacements in the given file"""
    renameRe = renamesSearchRe(replacements)
    outF = tempfile.NamedTemporaryFile(
        mode="w", suffix='', prefix='tmp', dir=os.path.dirname(filePath), delete=False, encoding='utf-8')
    didReplace = {}
    lineNr = 0
    with outF:
        with open(filePath, encoding='utf-8') as inF:
            while True:
                line = inF.readline()
                if not line:
                    break
                lineNr += 1
                ii = 0
                for m in renameRe.finditer(line):
                    old = m.group(1)
                    didReplace[old] = didReplace.get(old, []) + [lineNr]
                    outF.write(line[ii:m.start()])
                    outF.write(replacements[old])
                    ii = m.end()
                outF.write(line[ii:])
    if didReplace:
        bkPathBase = filePath + ".bk"
        bkPath = bkPathBase
        # non atomic (should really create if not there)
        ii = 0
        while os.path.exists(bkPath):
            ii += 1
            bkPath = bkPathBase + str(ii)
        os.rename(filePath, bkPath)
        os.rename(outF.name, filePath)
        print("%r: {" % filePath)
        for k,v in didReplace.items():
            print(k,'->',replacements[k],":",v)
        print("}\nBackup in ",bkPath)

if __name__ == "__main__":
    import sys
    renamesPath = os.path.join(os.path.dirname(__file__), "renames.txt")
    renames = loadRenamesFile(renamesPath)
    for f in sys.argv[1:]:
        try:
            replaceInFile(f, renames)
        except:
            logging.exception("handling file %s", f)
    print("DONE")
