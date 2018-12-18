import os, logging

from nomadcore.local_meta_info import InfoKindEl, InfoKindEnv, loadJsonFile

logger = logging.getLogger(__name__)
baseDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
metaInfoDir = os.path.normpath(os.path.join(baseDir, "../nomad_meta_info"))

allMetaInfo = {}
for dirpath, dirnames, filenames in os.walk(metaInfoDir):
    for filename in filenames:
        if not filename.endswith(".nomadmetainfo.json"):
            continue
        filepath = os.path.join(dirpath, filename)
        try:
            newEnv, warnings = loadJsonFile(filepath, dependencyLoader = None, extraArgsHandling = InfoKindEl.ADD_EXTRA_ARGS, uri = None)
            if warnings.get("hasWarnings", True):
                logger.warn("loading %s had the following warnings: %s" %
                            (filepath, warnings))
            allMetaInfo[filepath] = newEnv
        except:
            logger.exception("nomadcore.basic_meta_info could not load %s", filepath)
baseMetaInfo = allMetaInfo[os.path.join(metaInfoDir, "nomad_base.nomadmetainfo.json")]
metaMetaInfo = allMetaInfo[os.path.join(metaInfoDir, "meta_types.nomadmetainfo.json")]
