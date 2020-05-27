#!/bin/env python

import os
import tempfile
import re
import logging
from io import open
from io import StringIO


def loadRenamesFile(renamesPath, reverse=False):
    """Loads the content of a rename file"""
    renamesRe = re.compile(
        r"\s*(?P<old>[-a-zA-Z0-9_.%]+)\s*->\s*(?P<new>[a-zA-Z0-9_.]+)\s*$")
    renames = {}
    with open(renamesPath, encoding='utf-8') as f:
        while True:
            line = f.readline()
            if not line:
                break
            m = renamesRe.match(line)
            if m:
                if reverse:
                    renames[m.group("new")] = m.group("old")
                else:
                    renames[m.group("old")] = m.group("new")
            elif not line.isspace():
                logging.warn("Unexpected line %r in %s", line, renamesPath)
    return renames


def renamesSearchRe(renames):
    """creates a regular expression that matches the words that should be renamed"""
    res = StringIO()
    res.write(r"\b(")
    first = True
    for k in reversed(sorted(renames.keys())):
        if not first:
            res.write("|")
        else:
            first = False
        res.write(k.replace(".","\\."))
    res.write(r")(__|\b)")
    renameReStr = res.getvalue()
    return re.compile(renameReStr)


def replaceInFile(filePath, replacements):
    """performs the replacements in the given file"""
    renameRe = renamesSearchRe(replacements)
    outF = tempfile.NamedTemporaryFile(
        mode="w", suffix='', prefix='tmp', dir=os.path.dirname(filePath),
        delete=False, encoding='utf-8')
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
                    if m.group(2):
                        outF.write(m.group(2))
                    ii = m.end()
                outF.write(line[ii:])
    if didReplace:
        import datetime
        t=datetime.date.today()
        bkPathBase = "%s.%d-%02d-%02d" % (filePath, t.year, t.month, t.day)
        bkPath = bkPathBase + ".bk"
        # non atomic (should really create if not there)
        ii = 0
        while os.path.exists(bkPath):
            ii += 1
            bkPath = "%s-%d.bk" % (bkPathBase,ii)
        os.rename(filePath, bkPath)
        os.rename(outF.name, filePath)
        print("%r: {" % filePath)
        for k, v in didReplace.items():
            print(k, '->', replacements[k], ":", v)
        print("}\nBackup in ", bkPath)
    else:
        os.remove(outF.name)

if __name__ == "__main__":
    import sys,argparse
    renamesPath = os.path.join(os.path.dirname(__file__), "renames.txt")
    parser = argparse.ArgumentParser(description='Make replacements in files')
    parser.add_argument('--reverse', action='store_true',
                        help='Performs replacements in the reverse direction')
    parser.add_argument('--renames-file', nargs=1,
                        default=[renamesPath],
                        help='file containing the replacements to perform (from -> to)')
    parser.add_argument('toRename', metavar='P', nargs='+',
                        help='path to a file to rename')
    args = parser.parse_args()
    if not args.renames_file:
        baseRenames = loadRenamesFile(renamesPath, args.reverse)
    else:
        baseRenames = loadRenamesFile(args.renames_file[0], args.reverse)
    for f in args.toRename:
        try:
            basename = os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0]
            specificRenames = os.path.join(os.path.dirname(renamesPath), basename + ".renames")
            if not args.renames_file and os.path.exists(specificRenames):
                renames = loadRenamesFile(specificRenames, arg.reverse)
                renames.update(baseRenames)
            else:
                renames = baseRenames
            replaceInFile(f, renames)
        except:
            logging.exception("handling file %s", f)
    print("DONE")
