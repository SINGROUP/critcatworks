#!/usr/bin/env python
"""Loads the given info in the specified db"""
import sys, os, os.path, logging
basePath = os.path.realpath(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if not basePath in sys.path:
    sys.path.append(basePath)
from nomadcore.local_meta_info import loadJsonFile
from nomadcore.model_meta_info import InfoKind
from nomadcore.model_base import createEngine, createDB, createSession
from nomadcore.basic_meta_info import metaMetaInfo
from nomadscripts.calculate_meta_info_overrides import loadPath, loadRef
import sqlite3

def addToDb(parsed, session):
    """loads a nomadmetainfo from the given file"""
    for filepath, env in parsed['envs'].items():
        try:
            InfoKind.ensureInfoKindEnv(env, session)
        except:
            logging.exception("Error adding %s", filepath)
    for filepath, overrides in parsed['overrides'].items():
        pass # to do

if __name__ == "__main__":
    defaultPath = os.path.normpath(os.path.join(basePath, "../nomad_meta_info"))
    defaultDb = "sqlite:///" + os.path.normpath(os.path.join(basePath, "../default.sqlite"))

    usage = """{programName} [--db <dbString>] [--[no-]create-db] [--echo-sql] [--ref <repoRef> | --tags-matching <tagRegexp>]
    [--[no-]overrides] [--set-as-default] [--repo-path path] [path1 ...]

    loads the kindinfo.json files from basePath (defaults to {defaultPath})
    into the database reached by the given dbString (defaults to {defaultDb}).
    If repoRef is given then basePath has to be in a git repository, and the files
    are the the repoRef version.
    If tagRegexp is given all tags matching the given regexp are added to the db.
    By default overrides are not added to the db, but one can add them with --overrides.
    with --set-as-default it is ensured that the given version is see and the latest official version.
    --create-db creates the db, and --echo-sql echoes all sql commands.
    For sqlite connection if the file does not exist, then --create-db is automatically implied.
    """.format(programName = os.path.basename(sys.argv[0]),
               defaultDb = defaultDb, defaultPath = defaultPath)

    db = defaultDb
    iArg = 0
    nArgs = len(sys.argv)
    pathToLoad = []
    echoSql = False
    createDb = None
    overrides = False
    makeDefault = False
    refName = None
    tagRegExp = None
    repoPath = None

    while iArg + 1 < nArgs:
        iArg += 1
        arg = sys.argv[iArg]
        if arg == "--db":
            iArg += 1
            if iArg < nArgs:
                db = sys.argv[iArg]
        elif arg == "--create-db":
            createDb = True
        elif arg == "--no-create-db":
            createDb = False
        elif arg == "--echo-sql":
            echoSql = True
        elif arg == "--help":
            sys.stdout.write(usage)
            sys.exit(0)
        elif arg == "--ref":
            if tagRegExp:
                sys.stderr.write("\nError: defined both --ref and --tags-matching\n\n")
                sys.sterr.write(usage)
                sys.exit(1)
            elif ref:
                sys.stderr.write("\nError: multiple definitions of --ref\n\n")
                sys.sterr.write(usage)
                sys.exit(2)
            iArg += 1
            if iArg < nArgs:
                refName = sys.argv[iArg]
            else:
                sys.stderr.write("\
                sys.stderr.wError: missing reference after --ref\n\n")
                sys.sterr.write(usage)
                sys.exit(2)
        elif arg == "--tags-matching":
            if ref:
                sys.stderr.write("\nError: defined both --ref and --tags-matching\n\n")
                sys.sterr.write(usage)
                sys.exit(3)
            elif tagRegExp:
                sys.stderr.write("\nError: multiple definitions of --tags-matching\n\n")
                sys.sterr.write(usage)
                sys.exit(4)
            iArg += 1
            if iArg < nArgs:
                refName = sys.argv[iArg]
            else:
                sys.stderr.write("\nError: missing reference after --tags-matching\n\n")
                sys.sterr.write(usage)
                sys.exit(5)
        elif arg == '--overrides':
            overrides = True
        elif arg == '--no-overrides':
            overrides = False
        elif arg == '--set-as-default':
            makeDefault = True
        elif arg == "--repo-path":
            iArg += 1
            if iArg < len(sys.argv):
                repoPath = sys.argv[iArg]
            else:
                sys.stderr.write("\nError: missing path after --repo-path\n\n")
                sys.stderr.write(usage)
                sys.exit(5)
        elif os.path.exists(arg) or refName or tagRegExp:
            pathToLoad.append(arg)
        else:
            sys.stderr.write("\nError: expected a valid path\n\n")
            sys.stderr.write(usage)
            sys.exit(6)
    if not pathToLoad:
        pathToLoad.append(defaultPath)
    if db.startswith('sqlite://') and createDb is None:
        dbPath = db[len('sqlite://'):]
        if not os.path.exists(dbPath) and os.path.exists(os.path.dirname(dbPath)):
            createDb = True
    engine = createEngine(db, echo = echoSql)
    if createDb:
        createDB(engine)
    session = createSession(engine)
    if refName or tagRegExp:
        requireRepo = True
    else:
        requireRepo = False
    InfoKind.ensureInfoKindEnv(metaMetaInfo, session)
    if requireRepo and repoPath is None:
        repoPath = os.path.normpath(os.path.abspath(pathToLoad[0]))
        while not os.path.exists(os.path.join(repoPath, ".git")):
            if repoPath == '/':
                sys.stderr.write("Error: path {0!r} from is not within a git repo.\nCould not find .git directory, for bare repositories use the --repo-path option\n\n".format(pathToLoad[0]))
                sys.stderr.write(usage)
                sys.exit(10)
            repoPath = os.path.split(repoPath)[0]
    if refName:
        rev = repo.rev_parse(refName)
        revSha = oldRev.hexsha[:10]
        for path in pathToLoad:
            rPath = os.path.relpath(path, repoPath)
            try:
                env = loadRef(rev, rPath, loadOverrides = overrides)
                addToDb(env, session, )
            except:
                logging.exception("for path %s in revision %s (%s)", rPath, refName, revSha)
    elif tagRegExp:
        tagRe = re.compile(tagRegExp)
        repo = git.Repo(repoPath)
        for tag in repo.tags:
            if tagRe.match(tag.name):
                try:
                    ref = tag.ref
                except:
                    logging.exception("failed to resolve tag %s", tag.name)
                else:
                    for path in pathToLoad:
                        rPath = os.path.relpath(path, repoPath)
                        try:
                            env = loadRef(ref, rPath, loadOverrides = overrides)
                            addToDb(env, session)
                        except:
                            logging.exception("for path %s in tag %s (%s)", rPath, tag.name, ref.hexsha[:10])
    else:
        for path in pathToLoad:
            env = loadPath(path, loadOverrides = overrides)
            addToDb(env, session)
    session.commit()
