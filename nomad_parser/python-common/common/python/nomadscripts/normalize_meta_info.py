#!/usr/bin/env python
"""
Normalizes i.e. checks consistency, completness, absence of duplicates, and reformats meta_info data
"""
from builtins import str
import sys, os, os.path, datetime
basePath = os.path.realpath(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if not basePath in sys.path:
    sys.path.append(basePath)
from nomadcore.local_meta_info import InfoKindEnv, loadJsonFile, InfoKindEl
from nomadcore.json_support import jsonIndentD, jsonCompactD
from nomadcore.compact_sha import sha512

def normalizeFile(path, checkOnly, force, addGid, addSubGids, extraArgsHandling, s = sys.stderr.write):
    """Normalizes the meta info file at the given path"""
    env, warnings = loadJsonFile(path, extraArgsHandling = extraArgsHandling)
    hasErrors = False
    ovr = warnings.get("overrides")
    relpath = os.path.relpath(path)
    if ovr:
        hasErrors = True
        s("Error: loading {0!r}, found duplicates:\n".format(relpath))
        for old, new in ovr:
            jsonIndentD(old, s, extraIndent = 4)
            s(" => ")
            jsonIndentD(new, s, extraIndent = 4)
            s("\n")
    dup = warnings.get("duplicates")
    if dup:
        hasErrors = True
        s("Error: loading {0!r},\n    the depenendencies already define kinds with the following names:\n    ".format(relpath))
        s(str(dup))
        s("\n")
    hid = warnings.get("hidden")
    if hid:
        hasErrors = True
        s("Error: loading {0!r},\n    the depenendencies already define different kinds with the same names:\n    ".format(
                relpath))
        s(str([x.toDict() for x in hid]))
        s("\n")
    res = {"path":path, "env": env, "hasErrors": hasErrors, "hasChanges": False}
    try:
        env1 = InfoKindEnv(name = env.name, description = env.description, infoKinds=list(env.infoKinds.values()), path=env.path, uri=env.uri, deps=env.deps)
        env1.calcGids()
        validGids = True
    except:
        validGids = False
        hasErrors = True
        s("Error: parsing {0!r} failed to re-generate gids: {1}\n".format(relpath, sys.exc_info()[1]))
    if (hasErrors and not force):
        return res
    if env.gids and validGids:
        changes = []
        noDepNames = env.noDepNames()
        for k, v in env.gids.items():
            if k in noDepNames and k in env1.gids and env1.gids.get(k) != v:
                changes.append(k)
        if changes:
            s("Info: {0!r} regenerating gids changed for: ".format(relpath))
            jsonCompactD(changes, s)
            s("\n")
    if validGids:
        env = env1
    oldSha = None
    with open(path) as f:
        sha = sha512()
        while True:
            block = f.read(8*1024)
            if not block:
                break
            sha.update(block)
        oldSha = sha.b64digest()
    newS = sha512()
    try:
        env.serialize(newS.update, subGids = addSubGids, selfGid = addGid)
        newSha = newS.b64digest()
    except:
        hasErrors = True
        res["hasErrors"] = hasErrors
        s("Error: serialization of {0!r} failed: {1}".format(relpath, sys.exc_info()[1]))
        return res
    backup = None
    if (oldSha == newSha):
        s(repr(path))
        s(": OK, no change\n")
    else:
        backupBase = path + "-" + datetime.date.today().isoformat()
        extra = ""
        n = 0
        while os.path.exists(backupBase + extra + ".bk"):
            n += 1
            extra = "-" + str(n)
        backup = backupBase + extra + ".bk"
        if not checkOnly:
            os.rename(path, backup)
            with open(path, "wb") as f:
                env.serialize(f.write, subGids = addSubGids, selfGid = addGid)
        s(repr(path))
        s(": %s, has changes, backup in " % ("OK" if not hasErrors else "tying to fix ERRORS"))
        s(os.path.basename(backup))
        s("\n")
        res["hasChanges"] = True
        res["backup"] = backup
    return res

def normalizePath(path, checkOnly, force, addGid, addSubGids, extraArgsHandling, s = sys.stderr.write):
    """Normalizes all the nomadmetainfo.json files at the given path"""
    if (os.path.isfile(path)):
        return normalizeFile(path, checkOnly=checkOnly, force=force, addGid=addGid, addSubGids=addSubGids, extraArgsHandling=extraArgsHandling, s=s)
    envNames = {}
    allNames = set()
    hasErrors = False
    hasChanges = False
    clashes = set()
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            if not filename.endswith(".nomadmetainfo.json"):
                if not filename.endswith(".bk") and not filename.startswith("."):
                    s("skipping ")
                    s(filepath)
                    s("\n")
                continue
            fileRes = normalizeFile(filepath, checkOnly=checkOnly, force=force,
                    addGid=addGid, addSubGids=addSubGids, extraArgsHandling=extraArgsHandling, s=s)
            hasChanges = hasChanges or fileRes.get("hasChanges", False)
            hasErrors = hasErrors or fileRes.get("hasErrors", False)
            newNames = fileRes.get("env").noDepNames()
            envNames[filepath] = newNames
            newClashes = allNames.intersection(newNames)
            allNames.update(newNames)
            clashes.update(newClashes)
    detailedClashes = {}
    for filePath, names in envNames.items():
        currentClashes = names.intersection(clashes)
        if currentClashes:
            detailedClashes[filePath]=currentClashes
    if clashes:
        hasErrors = True
        s(repr(path))
        s(", detected clashes:\n")
        for filePath, fileClashes in detailedClashes.items():
            s("    ")
            s(filePath)
            s(":\n        ")
            jsonCompactD(list(fileClashes), s)
            s("\n")
    return {"path": path, "hasErrors": hasErrors, "hasChanges": hasChanges}

def normalizePaths(paths, checkOnly, force, addGid, addSubGids, extraArgsHandling, s = sys.stderr.write):
    """Normalizes all the given paths"""
    hasErrors = False
    hasChanges = False
    messages = []
    for path in paths:
        messages.append(path)
        res = normalizePath(path, checkOnly=checkOnly, force=force,
                addGid=addGid, addSubGids=addSubGids, extraArgsHandling=extraArgsHandling, s=s)
        errors = res.get("hasErrors", True)
        hasErrors = hasErrors or errors
        if not errors:
            messages.append(": OK")
        else:
            messages.append(": ERROR")
        changes = res.get("hasChanges", False)
        hasChanges = hasChanges or changes
        if changes:
            messages.append(", with changes\n")
        else:
            messages.append(", no changes\n")
    s("\n")
    for m in messages:
        s(m)
    s("\n")
    if not hasErrors:
        s("OK")
    else:
        s("ERROR")
    if hasChanges:
        s(", with changes\n")
        if checkOnly:
            s("    (run performed with --check-only, no changes were committed to the filesystem)\n")
    else:
        s(", no changes\n")
    return {"paths": paths, "hasErrors": hasErrors, "hasChanges": hasChanges}

if __name__ == "__main__":
    defaultDir = os.path.normpath(os.path.join(basePath, "../../../nomad-meta-info/meta_info/nomad_meta_info"))
    usage = """usage: {command} [--check-only] [--force] [--add-gid] [--add-sub-gids] [--keep-extra-args] [--remove-extra-args] [--error-if-extra-args] [file/or/dir ...]

    normalizes the InfoKinds file/or/dir (that defalts to {dir})
""".format(command=os.path.basename(sys.argv[0]), dir=defaultDir)

    checkOnly = False
    force = False
    addGid = False
    addSubGids = False
    paths = []
    extraArgsHandling = InfoKindEl.ADD_EXTRA_ARGS

    for arg in sys.argv[1:]:
        if arg == "--check-only":
            checkOnly = True
        elif arg == "--force":
            force = True
        elif arg == "--add-gid":
            addGid = True
        elif arg == "--add-sub-gids":
            addSubGids = True
        elif arg == "--add-extra-args":
            extraArgsHandling = InfoKindEl.ADD_EXTRA_ARGS
        elif arg == "--remove-extra-args":
            extraArgsHandling = InfoKindEl.IGNORE_EXTRA_ARGS
        elif arg == "--error-if-extra-args":
            extraArgsHandling = InfoKindEl.RAISE_IF_EXTRA_ARGS
        elif arg == "--help":
            sys.stdout.write(usage)
            sys.exit(0)
        elif os.path.exists(arg):
            paths.append(arg)
        else:
            sys.stdout.write("Error: invalid argument ")
            sys.stdout.write(repr(arg))
            sys.stdout.write("\n")
            sys.stdout.write(usage)
            sys.exit(1)
    if not paths:
        paths.append(defaultDir)
    normalizePaths(paths, checkOnly=checkOnly, force=force, addGid=addGid, addSubGids=addSubGids, extraArgsHandling=extraArgsHandling, s = sys.stderr.write)
