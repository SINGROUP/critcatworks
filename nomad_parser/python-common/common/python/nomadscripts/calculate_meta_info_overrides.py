#!/usr/bin/env python
from __future__ import print_function
from past.builtins import cmp
from builtins import range
from builtins import object
import sys, os, os.path, datetime, logging, json
basePath = os.path.realpath(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if not basePath in sys.path:
    sys.path.append(basePath)
from nomadcore.local_meta_info import InfoKindEnv, InfoKindEl, loadJsonFile, loadJsonStream
from nomadcore.json_support import jsonIndentF, jsonCompactF
import git

class GitTreeDependencySolver(object):
    def __init__(self, tree, basePath, context):
        self.tree = tree
        self.basePath = basePath
        self.context = context
        self.deps = {}

    def __call__(self, infoKindEnv, source, dep):
        if "relativePath" not in dep:
            raise Exception('Invalid dependency for relativeDependencySolver there must be a relativePath')
        dPath = os.path.normpath(os.path.join(self.basePath, dep.get('relativePath')))
        if dPath in self.deps:
            return self.deps[dPath]
        depInfo = None
        depIKEnv = InfoKindEnv(name = os.path.basename(dPath), path = dPath, dependencyLoader=infoKindEnv.dependencyLoader)
        self.deps[dPath] = depIKEnv
        try:
            f=self.tree[dPath].data_stream
        except:
            print(self.tree.hexsha, ", t.path", self.tree.path, ", bpath:", self.basePath, ", dPath:", dPath)
            raise Exception("Failed to resolve dependency {0!r} in {context}, due to exception {1}".format(dep, sys.exc_info()[1], context = self.context))
        depInfo = json.load(f)
        if depInfo:
            depIKEnv.fromJsonList(depInfo, source = { 'path': dPath }, dependencyLoad = True)
        return depIKEnv

def loadPath(path, loadOverrides = False):
    """loads all nomadmetainfo.json files within the given path.
    If loadOverrides is true then also nomadmetainfooverrides.json files are loaded"""
    allEnvs = {}
    allOverrides = {}
    hasWarnings = False
    if os.path.isdir(path):
        for dirpath, dirnames, filenames in os.walk(path):
            for filename in filenames:
                if filename.endswith(".nomadmetainfo.json"):
                    filepath = os.path.join(dirpath, filename)
                    env, warnings = loadJsonFile(filepath, extraArgsHandling = InfoKindEl.ADD_EXTRA_ARGS)
                    if warnings.get("duplicates", False) or warnings.get("overrides", False):
                        hasWarnings = True
                        logging.warn("{0!r}: warnings in loadJsonFile:{1}\n".format(filepath, warnings))
                    allEnvs[os.path.relpath(filepath, path)] = env
                elif loadOverrides and filename.endswith(".nomadmetainfooverrides.json"):
                    try:
                        with open(filename) as f:
                            overrides = json.load(f)
                    except:
                        hasWarnings = True
                        logging.exception("Error loading %r", filename)
                    else:
                        allOverrides[os.path.relpath(filepath, path)] =  overrides
    return {"hasWarnings": hasWarnings, "envs": allEnvs, "overrides": allOverrides}

def loadRef(rev, path, loadOverrides = False):
    """loads all nomadmetainfo.json files within path from repo reference"""
    allEnvs = {}
    allOverrides = {}
    hasWarnings = False
    if path == ".":
        path = ""
    if path:
        tree = rev.tree[path]
    else:
        tree = rev.tree
    def loadBlobHasWarnings(obj):
        hasW = False
        if obj.type != 'blob':
            return False
        rPath = os.path.relpath(obj.path, path)
        if obj.name.endswith(".nomadmetainfo.json"):
            env, warnings = loadJsonStream(obj.data_stream, GitTreeDependencySolver(rev.tree, os.path.dirname(obj.path), obj.path),
                                           filePath = obj.path)
            allEnvs[rPath] = env
            dup = warnings.get("duplicates")
            hid = warnings.get("hidden")
            if dup or hid:
                hasW = True
                logging.warn("loading {path} of revision {rev}: {warnings}\n".format(
                    path = obj.path, rev = rev.name_rev, warnings = warnings))
        elif obj.name.endswith(".nomadmetainfooverrides.json"):
            try:
                overrides = json.load(obj.data_stream)
            except:
                hasW = True
                logging.exception("Error loading %r in revision %s", obj.path, rev.hexsha[:10])
            else:
                allOverrides[rPath] = overrides
        return hasW
    if tree.type == 'blob':
        if loadBlobHasWarnings(tree):
            hasWarnings = True
    else:
        for obj in tree.traverse(lambda i, d: i.type == 'blob' and (i.path.endswith(".nomadmetainfo.json") or
                                                                    (loadOverrides and i.path.endswith(".nomadmetainfooverrides.json")))):
            if loadBlobHasWarnings(obj):
                hasWarnings = True
    return {"hasWarnings": hasWarnings, "envs": allEnvs, "overrides": allOverrides}

def insertDicList(dict, key, value):
    if key in dict:
        dict[key].append(value)
    else:
        dict[key] = [value]

def prepareValues(values, name):
    "prepares the values computing gids"
    allNames = {}
    hasWarnings = False
    for path, env in values.items():
        try:
            env.calcGids()
        except:
            hasWarnings = True
            logging.warn("Error calculating gids of {path!r} {name}: {exc}\n".format(path = path, name = name, exc = sys.exc_info()[1]))
        noDepNames = env.noDepNames()
        for name, gid in env.gids.items():
            if name in noDepNames:
                insertDicList(allNames, name, (path, gid))
    return {"allNames": allNames, "envs": values, "hasWarnings": hasWarnings}

def cmpOverrides(o1,o2):
    c = cmp(o1["name"], o2["name"])
    if c != 0:
        return c
    return cmp(o1, o2)

def buildOverrides(oldValues, newValues):
    oldV = prepareValues(oldValues, "(old)")
    newV = prepareValues(newValues, "(new)")
    oldNames = oldV["allNames"]
    newNames = newV["allNames"]
    overrides = []
    complexOverrides = []
    additions = []
    removals = []
    hasWarnings = oldV["hasWarnings"] or newV["hasWarnings"]
    nMatched = 0
    for name, oldUses in oldNames.items():
        if not name in newNames:
            removals.append({"name": name, "oldGids": [x[1] for x in oldUses], "newGids": []})
            continue
        newUses = newNames[name]
        newNoMatch = list(newUses)
        oldNoMatch = []
        for oldPath, oldGid in oldUses:
            found = -1
            for i in range(len(newNoMatch)):
                newPath, newGid = newNoMatch[i]
                if newGid == oldGid:
                    del newNoMatch[i]
                    found = i
                    break
            if found != -1:
                nMatched += 1
            else:
                oldNoMatch.append((oldPath, oldGid))
        if not oldNoMatch and not newNoMatch:
            continue
        if len(oldNoMatch) == 1 and len(newNoMatch) == 1:
            overrides.append({"name": name, "oldGid": oldNoMatch[0][1], "newGid": newNoMatch[0][1]})
            continue
        # try path matching
        iOld = 0
        while iOld < len(oldNoMatch):
            oldPath, oldGid = oldNoMatch[iOld]
            found = -1
            for iNew in range(len(newNoMatch)):
                newPath, newGid = newNoMatch[iNew]
                if newPath == oldPath:
                    overrides.append({"name": name, "oldGid": oldGid, "newGid": newGid})
                    del newNoMatch[iNew]
                    del old
                    found = iNew
                    break
            if found != -1:
                del oldNoMatch[iOld]
            else:
                iOld += 1
        if len(oldNoMatch) == 1 and len(newNoMatch) == 1:
            overrides.append({"name": name, "oldGid": oldNoMatch[0][1], "newGid": newNoMatch[0][1]})
        if not oldNoMatch and not newNoMatch:
            continue
        elif oldNoMatch and not newNoMatch:
            removals.append({"name": name, "oldGids": [x[1] for x in oldNoMatch], "newGids": []})
        elif not oldNoMatch and newNoMatch:
            additions.append({"name": name, "oldGids": [], "newGids": [x[1] for x in newNoMatch]})
        else: # oldNoMatch and newNoMatch
            complexOverrides.append({"name": name, "oldGids": oldNoMatch, "newGids": newNoMatch[0][1]})
    for name in newNames.keys():
        if not name in oldNames:
            additions.append({"name": name, "oldGids": [], "newGids": [x[1] for x in newNames[name]]})
    overrides.sort(cmpOverrides)
    complexOverrides.sort(cmpOverrides)
    additions.sort(cmpOverrides)
    removals.sort(cmpOverrides)
    return {
        "oldNames": oldNames,
        "newNames": newNames,
        "overrides": overrides,
        "complexOverrides": complexOverrides,
        "additions": additions,
        "removals": removals,
        "hasWarnings": hasWarnings,
        "nMatched": nMatched
    }



if __name__ == "__main__":
    overridesDir = os.path.normpath(os.path.join(basePath,"../../../nomad-meta-info/meta_info/nomad_meta_info_overrides"))
    defaultPath = os.path.normpath(os.path.join(basePath,"../../../nomad-meta-info/meta_info/nomad_meta_info"))
    usage = """usage: {command} [--check-only] [--old-ref <ref1=HEAD>] [--old-path <path=None>] [--help]
    [--new-ref <ref2=None>] [--new-path <path=basePath>] [--repo-path <repoPath=None>]
    [--overrides-dir <overridesDir>] [--no-clobber] [--verbose] [<basePath>]

Generates the overrides for the kind info with the same name but different sha from the old version
to the new version.
Old version can be given as a revision in git (defaults to HEAD) with --old-ref,  or a path with
--old-path, likewise the new version can be given as a revision with --new-ref or new-path.
basePath (which defaults to {defaultPath}) gives the path of which the revisions should be compared (might be a
subdirectory or file of the repo) and if no --new-ref is given the version checked out there is
directly the new vesion.
The repository path by default is found looking from path upward for a .git repository, but can
also be given explicitly with --repo-path.
overridesDir defaults to {overridesDir}
By default output goes to stdout, but an explicit --out-file can also be given
""".format(command = os.path.basename(sys.argv[0]), overridesDir = overridesDir, defaultPath = defaultPath)
    oldRef = None
    oldPath = None
    newRef = None
    newPath = None
    repoPath = None
    path = None
    checkOnly = False
    noClobber = False
    verbose = False

    iArg = 1
    while iArg < len(sys.argv):
        arg = sys.argv[iArg]
        if arg == "--help":
            sys.stderr.write(usage)
            sys.exit(0)
        elif arg == "--check-only":
            checkOnly = True
        elif arg == "--no-clobber":
            noClobber = True
        elif arg == "--verbose":
            verbose = True
        elif arg == "--old-ref":
            iArg += 1
            if iArg < len(sys.argv):
                oldRef = sys.argv[iArg]
            else:
                sys.stderr.write("Error: missing reference after --old-ref\n\n")
                sys.stderr.write(usage)
                sys.exit(1)
        elif arg == "--old-path":
            iArg += 1
            if iArg < len(sys.argv):
                oldPath = sys.argv[iArg]
            else:
                sys.stderr.write("Error: missing path after --old-path\n\n")
                sys.stderr.write(usage)
                sys.exit(2)
        elif arg == "--new-ref":
            iArg += 1
            if iArg < len(sys.argv):
                newRef = sys.argv[iArg]
            else:
                sys.stderr.write("Error: missing reference after --new-ref\n\n")
                sys.stderr.write(usage)
                sys.exit(3)
        elif arg == "--new-path":
            iArg += 1
            if iArg < len(sys.argv):
                newPath = sys.argv[iArg]
            else:
                sys.stderr.write("Error: missing path after --new-path\n\n")
                sys.stderr.write(usage)
                sys.exit(4)
        elif arg == "--repo-path":
            iArg += 1
            if iArg < len(sys.argv):
                repoPath = sys.argv[iArg]
            else:
                sys.stderr.write("Error: missing path after --repo-path\n\n")
                sys.stderr.write(usage)
                sys.exit(5)
        elif arg == "--overrides-dir":
            iArg += 1
            if iArg < len(sys.argv):
                overridesDir = sys.argv[iArg]
            else:
                sys.stderr.write("Error: missing reference after --overrides-dir\n\n")
                sys.stderr.write(usage)
                sys.exit(7)
        elif path is None:
            path = arg
        else:
            sys.stderr.write("Error: unexpected extra argument {0!r}\n\n".format(arg))
            sys.stderr.write(usage)
            sys.exit(8)
        iArg += 1
    if path is None:
        path = defaultPath
    requiresRepo = True
    if oldPath is None and oldRef is None:
        oldRef = "HEAD"
    if not newRef and not newPath:
        newPath = path
    if oldPath and not newRef:
        requiresRepo = False
        if oldRef:
            sys.stderr.write("Error: both --old-path and --old-ref were given this is not supported\n\n")
            sys.stderr.write(usage)
            sys.exit(9)
    if requiresRepo and repoPath is None:
        repoPath = os.path.normpath(os.path.abspath(path))
        while not os.path.exists(os.path.join(repoPath, ".git")):
            if repoPath == '/':
                sys.stderr.write("Error: path {0!r} from is not within a git repo.\nCould not find .git directory, for bare repositories use the --repo-path option\n\n".format(arg))
                sys.stderr.write(usage)
                sys.exit(10)
            repoPath = os.path.split(repoPath)[0]
    if repoPath:
        repo = git.Repo(repoPath)
        if os.path.abspath(path).startswith(repoPath) and os.path.exists(os.path.abspath(path)):
            relPath = os.path.relpath(os.path.normpath(os.path.abspath(path)), repoPath)
        else:
            relPath = path
        if requiresRepo and relPath.startswith("../"):
            sys.stderr.write("Error: basePath should be in the repo, i.e. within {0!r}\n\n".format(repoPath))
            sys.stderr.write(usage)
            sys.exit(11)
    oldName = ""
    if oldRef:
        oldRev = repo.rev_parse(oldRef)
        oldName = oldRev.hexsha[:10]
        oldValues = loadRef(oldRev, relPath)
        sys.stdout.write("<old>: {relPath} in {oldRef}({oldName}) of repo at {repoPath}\n".format(
            relPath = relPath, oldRef = oldRef, oldName = oldName, repoPath = repoPath))
    else:
        sys.stdout.write("<old>: {oldPath}\n".format(oldPath = oldPath))
        oldValues = loadPath(oldPath)
    newName = ""
    if newRef:
        newRev = repo.rev_parse(newRef)
        newName = newRev.hexsha[:10]
        newValues = loadRef(newRev, relPath)
        sys.stdout.write("<new>: {relPath} in {newRef}({newName}) of repo at {repoPath}\n".format(
            relPath = relPath, newRef = newRef, newName = newName, repoPath = repoPath))
    else:
        sys.stdout.write("<new>: {newPath}\n".format(newPath = newPath))
        newValues = loadPath(newPath)
    res = buildOverrides(oldValues["envs"], newValues["envs"])
    if verbose:
        jsonIndentF(res, sys.stdout, check_circular = False)
        sys.stdout.write("\n")
    sys.stdout.flush()
    overrides = res["overrides"]
    if overrides and not checkOnly:
        overridesName = os.path.join(overridesDir, "{old}_{new}_{date}.nomadmetainfooverrides.json"
                .format(old = oldName, new = newName, date = datetime.date.today().isoformat()))
        n = 1
        while noClobber and os.path.exists(overridesName):
            overridesName = os.path.join(overridesDir, "{old}_{new}_{date}_{n}.nomadmetainfooverrides.json"
                    .format(old = oldName, new = newName, date = datetime.date.today().isoformat()), n = n)
            n += 1
        if os.path.exists(overridesName):
            overwriteStr = "Overwriting"
        else:
            overwriteStr = "Writing"
        with open(overridesName, "w") as f: # use os.open to guarantee no clashes seems like an overkill
            jsonIndentF(overrides, f, check_circular = False)
        sys.stdout.write("{0} overrides in {1!r}\n".format(overwriteStr, os.path.abspath(overridesName)))
    hasWarnings = oldValues["hasWarnings"] or newValues["hasWarnings"] or res["hasWarnings"]
    removals = res["removals"]
    if removals:
        sys.stdout.write("removals: ")
        jsonCompactF([x["name"] for x in removals], sys.stdout)
        sys.stdout.write("\n")
    additions = res["additions"]
    if additions:
        sys.stdout.write("additions: ")
        jsonCompactF([x["name"] for x in additions], sys.stdout)
        sys.stdout.write("\n")
    complexOverrides = res["complexOverrides"]
    if complexOverrides:
        hasWarnings = True
        sys.stdout.write("Warning: complexOverrides: ")
        jsonCompactF([x["name"] for x in complexOverrides], sys.stdout)
        sys.stdout.write("\n")
    sys.stdout.write("#overrides:{nov:3} #matches:{nm:3} #complexOverrides:{ncov:3} #additions:{nadd:3} #removals:{nrm}\n"
            .format(nm = res["nMatched"], nov = len(overrides), ncov = len(complexOverrides), nadd = len(additions), nrm = len(removals)))
    if not hasWarnings:
        sys.exit(0)
    else:
        sys.stderr("Had warnings\n")
        sys.exit(1)
