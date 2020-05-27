"""
This module contains base classes for building parsers for the NoMaD project
with an object-oriented approach.
"""
from builtins import str
from builtins import object

import os
import copy
import numpy as np
import logging
from abc import ABCMeta, abstractmethod
from nomadcore.unit_conversion import unit_conversion
from nomadcore.simple_parser import mainFunction
from nomadcore.local_backend import LocalBackend
from nomadcore.local_meta_info import load_metainfo
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import extractOnCloseTriggers, extractOnOpenTriggers
from nomadcore.caching_backend import ActiveBackend
import nomadcore.ActivateLogging
from future.utils import with_metaclass
logger = logging.getLogger("nomad")


class ParserInterface(with_metaclass(ABCMeta, object)):
    """This class provides an interface for parsing. The end-user will
    typically only interact with an object that derives from this class. All
    the input is given to this class and the parsing is done by calling the
    parse() method. The parsing output is determined by the backend object that
    is given in the constructor as a dependency. If no backend is specified, a
    local backend that outputs results into a dictionary will be used.

    This class controls if a specific parser version should be instantiated and
    provides basic information about the parser, such as the parser-info
    -dictionary and the name of the metainfo file.

    Attributes:
        main_parser: Object that actually does the parsing and is
            setup by this class based on the given contents.
        parser_context: A wrapper class for all the parser related information.
            This is contructed here and then passed onto the different
            subparsers.
    """
    metainfo_env = None

    def __init__(self, metainfo_to_keep=None, backend=None, default_units=None, metainfo_units=None, debug=False, log_level=logging.ERROR, store=True):
        """
    Args:
        main_file: A special file that can be considered the main file of the
            calculation.
        metainfo_to_keep: A list of metainfo names. This list is used to
            optimize the parsing process as optimally only the information
            relevant to these metainfos will be parsed.
        backend: An object to which the parser will give all the parsed data.
            The backend will then determine where and when to output that data.
        """
        self.debug = debug
        logger.setLevel(log_level)
        self.store = store
        self.initialize(metainfo_to_keep, backend, default_units, metainfo_units)

    def initialize(self, metainfo_to_keep, backend, default_units, metainfo_units):
        """Initialize the parser with the given environment.
        """
        self.parser_context = ParserContext()
        self.parser_context.metainfo_to_keep = metainfo_to_keep
        self.parser_context.file_service = FileService()
        self.parser_context.cache_service = CacheService()
        self.parser_context.parser_info = self.get_parser_info()
        self.main_parser = None

        # Setup the metainfo environment. All parsers that inherit from this
        # class will have a static class attribute that will store the metainfo
        # environment. This way every instance of a parser doesn't have to load
        # the environment separately because it is identical for each instance.
        if type(self).metainfo_env is None:
            metainfo_env, warn = load_metainfo(self.get_metainfo_filename())
            type(self).metainfo_env = metainfo_env
            self.parser_context.metainfo_env = metainfo_env
        else:
            self.parser_context.metainfo_env = type(self).metainfo_env

        # Initialize the backend. Use local backend if none given
        if backend is not None:
            self.parser_context.super_backend = backend(type(self).metainfo_env)
        else:
            self.parser_context.super_backend = LocalBackend(type(self).metainfo_env, debug=self.debug, store=self.store)

        # Check the list of default units
        default_unit_map = {}
        if default_units is not None:
            for unit in default_units:
                dimension = unit_conversion.ureg(unit).dimensionality
                old_value = default_unit_map.get(str(dimension))
                if old_value is not None:
                    raise LookupError("You can only specify one default value per dimension in the 'default_units' list. There are two different units given for the dimension '{}'".format(dimension))
                default_unit_map[str(dimension)] = unit

        # Check the list of metainfo units
        if metainfo_units is not None:
            for metaname, unit in metainfo_units.items():

                # Check that the unit is OK
                unit_conversion.ureg(unit)

                # Check that the metaname is OK
                meta = ParserInterface.metainfo_env.infoKinds.get(metaname)
                if meta is None:
                    raise KeyError("The metainfo name '{}' could not be found. Check for typos or try updating the metainfo repository.".format(metaname))

        # Save the default units
        self.parser_context.default_units = default_unit_map
        self.parser_context.metainfo_units = metainfo_units

    @abstractmethod
    def setup_version(self):
        """Deduce the version of the software that was used and setup a correct
        main parser. The main parser should subclass MainParser and be stored
        to the 'main_parser' attribute of this class. You can give the
        'parser_context' wrapper object in the main parser constructor to pass
        all the relevant data for it.
        """
        pass

    # @abstractmethod  # Keep this way for a file in order to avoid breaking
    def setup_main_parser(self, version_dictionary):
        """Setup the attribute 'main_parser' which contains the object that
        starts the parsing from the main file. The main parser should be
        decided based on the information contained on the given dictionary.

        Args:
            version_dictionary: Dictionary containing version information. It
            can include the version number, run type, etc.
        """
        pass

    def get_mainfile_regex(self):
        """Used to return the regular expression that is used to match the main
        file.

        Returns:
            str: regular expression as a string.
        """

    @abstractmethod
    def get_metainfo_filename(self):
        """This function should return the name of the metainfo file that is
        specific for this parser. When the parser is started, the metainfo
        environment is loaded from this file that is located within a separate
        repository (nomad-meta-info).

        Returns:
            A string containing the metainfo filename for this parser.
        """
        return None

    @abstractmethod
    def get_parser_info(self):
        """This function should return a dictionary containing the parser info.
        This info is printed to the JSON backend. it should be of the form:

            {'name': 'softwarename-parser', 'version': '1.0'}

        Returns:
            A dictionary containing information about this parser.
        """
        return None

    def parse(self, main_file):
        """Starts the actual parsing process, and outputs the results to the
        backend specified in the constructor.
        """
        # Check that the main file exists
        if not os.path.isfile(main_file):
            raise ValueError(
                "Couldn't find the main file '{}'. Check that the path is valid "
                "and the file exists on this path.".format(main_file)
            )

        self.parser_context.main_file = main_file
        self.setup_version()
        if not self.main_parser:
            logger.error("The main parser has not been set up.")

        self.main_parser.parse(main_file)

        # If using a local backend, the results will have been saved to a
        # separate results dictionary which should be returned.
        try:
            return self.parser_context.super_backend.results
        except AttributeError:
            return None


class FileService(object):
    """Provides the interface to accessing files related to a calculation.

    Before you can use the service you must setup the root path, where every
    file related to this calculation resides. All queries filepaths will be
    handled as relative to this root folder.

    You can also setup ids that point to a certain path. This helps in querying
    the files as you don't have to remember the exact path and you can store
    paths for later use.

    Used to map file paths to certain ID's. This helps in setting up the
    Secondary parsers as you can associate file paths to simpler ID's that are
    easier to use.

    Attributes:
        root_folder: A path to the root folder
        _file_ids: A dictionary containing the mapping between file ids and filepaths
        _file_handles: A "private" dictionary containing the cached file handles
        _file_contents: A "private" dictionary containing the cached file contents
        _file_sizes: A "private" dictionary containing the cached file sizes
    """
    def __init__(self, root_folder=None):
        """
        Args:
            root_folder: a path to the root folder as a string.
        """
        self.map_id_to_path = {}
        if root_folder is not None:
            self.setup_root_folder(root_folder)

    def setup_root_folder(self, root_folder):
        """Setup the path to the root folder. Every filepath you set or get
        through this service should be relative to this root folder.
        """
        if os.path.isdir(root_folder):
            self.root_folder = root_folder
        else:
            raise IOError("Could not find the folder '{}'".format(root_folder))

    def get_absolute_path_to_file(self, relative_path):
        """
        Returns:
            If the given file is found, returns the absolute path to it. Return
            none if no file with the given path can be found.
        """
        path = os.path.join(self.root_folder, relative_path)
        if os.path.isfile(path):
            return path
        else:
            return None

    def get_file_by_id(self, file_id):
        """
        Returns:
            Handle to the file. Return none if no file with the given id has
            been set.
        """
        path = self.map_id_to_path.get(file_id)
        if path is None:
            logger.error("The id '{}' has no path associated with it.".format(file_id))
            return None
        else:
            return path

    def set_file_id(self, path, file_id):
        """Used to map a simple identifier string to a file path. When a file
        id has been setup, you can easily access the filepath anywhere in the
        code (FileService is shared via parser_context) by calling
        get_file_by_id().
        """
        if path is None:
            return None

        # Check if there is an old definition
        old = self.map_id_to_path.get(file_id)
        if old is not None:
            raise LookupError("The path '{}' is already associated with id '{}'".format(old, file_id))
        # Check that the file exists
        path = os.path.join(self.root_folder, path)
        if not os.path.isfile(path):
            logger.warning("Could not set the id for file '{}' as it cannot be found.".format(path))
            return None
        else:
            self.map_id_to_path[file_id] = path
            return path


class RegexService(object):
    """
    Stores basic regex definitions that can be reused on multiple parsers.
    """
    def __init__(self):
        self.float = "-?\d+\.\d+(?:[ED](?:\+|-)\d+)?"  # Regex for a floating point value
        self.int = "-?\d+"  # Regex for an integer
        self.word = "[\S]+"  # Regex for a single word. Can contain anything else but whitespace
        self.letter = "[^\W\d_]"  # Regex for a single alphabetical letter
        self.eol = "[^\n]+"  # Regex for a single alphabetical letter


class AbstractBaseParser(with_metaclass(ABCMeta, object)):
    """A base class for all objects that parse contents from files.

    When initialized with the parser_context, you can find the caching backend
    and the super backend as parameters of this object. These can then be used
    by the object to output parsing results.

    Attributes:
        file_path: Path to a file that is parsed by this class.
        parser_context: The ParserContext object that contains various in-depth
            information about the parsing environment.
        backend: The backend that will cache things according to the rules
            given in the main parser.
        super_backend: The final backend where things are forwarded to by the
            caching backend.
    """

    def __init__(self, parser_context):
        self.parser_context = parser_context
        self.backend = parser_context.caching_backend
        self.super_backend = parser_context.super_backend
        self.file_service = parser_context.file_service
        self.cache_service = parser_context.cache_service
        self.caching_levels = {}
        self.default_data_caching_level = CachingLevel.ForwardAndCache
        self.default_section_caching_level = CachingLevel.Forward
        self.on_close = {}
        self.on_open = {}

    def prepare(self):
        """This function will prepare everything for parsing.

        The onClose and onOpen callbacks are gathered and the ActiveBackend is
        initialized. You should call this function, or prepare things manually
        before trying to use push values to the backend.
        """
        # Gather the onClose and onOpen triggers
        for attr, callback in extractOnCloseTriggers(self).items():
            oldCallbacks = self.on_close.get(attr, None)
            if oldCallbacks:
                oldCallbacks.append(callback)
            else:
                self.on_close[attr] = [callback]
        for attr, callback in extractOnOpenTriggers(self).items():
            oldCallbacks = self.on_open.get(attr, None)
            if oldCallbacks:
                oldCallbacks.append(callback)
            else:
                self.on_open[attr] = [callback]

        # Initialize the Caching backend
        self.backend = ActiveBackend.activeBackend(
            metaInfoEnv=self.parser_context.metainfo_env,
            cachingLevelForMetaName=self.caching_levels,
            defaultDataCachingLevel=self.default_data_caching_level,
            defaultSectionCachingLevel=self.default_section_caching_level,
            onClose=self.on_close,
            onOpen=self.on_open,
            superBackend=self.parser_context.super_backend,
            default_units=self.parser_context.default_units,
            metainfo_units=self.parser_context.metainfo_units)

    def print_json_header(self):
        self.super_backend.fileOut.write("[")
        uri = "file://" + self.parser_context.main_file
        self.backend.startedParsingSession(uri, self.parser_context.parser_info)

    def print_json_footer(self):
        self.backend.finishedParsingSession("ParseSuccess", None)
        self.super_backend.fileOut.write("]\n")

    @abstractmethod
    def parse(self, filepath):
        """Used to do the actual parsing. Inside this function you should push
        the parsing results into the Caching backend, or directly to the
        superBackend. You will also have to open new sections, but remember
        that certain sections may already be opened by other parsers.
        """


class MainHierarchicalParser(AbstractBaseParser):
    """A base class for all parsers that parse a file using a hierarchy of
    SimpleMatcher objects.

    Attributes:
        root_matcher: The SimpleMatcher object at the root of this parser.
            caching_levels: A dicionary containing the caching options
            that the ActiveBackend will use for individual metanames.

            Example:
                self.caching_levels = {
                    'section_XC_functionals': CachingLevel.ForwardAndCache,
                }
        caching_levels: Dictionary that stores the caching levels for
            individual metanames.
        default_data_caching_level: A default caching level for data, i.e.
            metainfo with kindStr=type_document_content or type_dimension
        default_section_caching_level: A default caching level for sections.
        onClose: a dictionary of onClose triggers. These are added to the
            triggers that are already present in the parser context object.
        cm: An optional CommonMatcher object that is used to store common
            onClose triggers, SimpleMatchers, caching levels, etc.

    """
    def __init__(self, parser_context):
        """
        Args:
            file_path: Path to the main file as a string.
            parser_context: The ParserContext object that contains various
                in-depth information about the parsing environment.
        """
        super(MainHierarchicalParser, self).__init__(parser_context)
        self.root_matcher = None
        self.regexs = RegexService()
        self.cm = None
        self.super_context = self

    def parse(self, filepath):
        """Starts the parsing. By default uses the SimpleParser scheme, if you
        want to use something else or customize the process just override this
        method in the subclass.
        """
        mainFunction(
                mainFileDescription=self.root_matcher,
                metaInfoEnv=self.parser_context.metainfo_env,
                parserInfo=self.parser_context.parser_info,
                outF=self.parser_context.super_backend.fileOut,
                cachingLevelForMetaName=self.caching_levels,
                superContext=self.super_context,
                onClose=self.on_close,
                onOpen=self.on_open,
                default_units=self.parser_context.default_units,
                metainfo_units=self.parser_context.metainfo_units,
                superBackend=self.parser_context.super_backend,
                metaInfoToKeep=self.parser_context.metainfo_to_keep,
                mainFile=filepath)

    def startedParsing(self, fInName, parser):
        """Function is called when the parsing starts.

        Args:
            fInName: The name of the main file.
            parser: The SimpleParser object that has compiled the regex
                definitions inside the root_matcher attribute of this class.
        """
        self.parser_context.caching_backend = parser.backend
        self.parser_context.cache_service.parser_context = self.parser_context
        self.parser_context.super_backend = parser.backend.superBackend
        self.backend = self.parser_context.caching_backend
        self.super_backend = self.parser_context.super_backend
        if self.cm is not None:
            self.cm.backend = parser.backend

    def setup_common_matcher(self, common_matcher):
        """Used to setup the CommonMatcher object. This object will contain
        SimpleMatchers, onClose functions, caching levels, etc. that are common
        for many different HierarchicalParsers.

        Args:
            common_matcher: A CommonMatcher object from which to extract stuff.
        """
        self.cm = common_matcher
        self.on_close.update(common_matcher.getOnCloseTriggers())
        self.on_open.update(common_matcher.getOnOpenTriggers())
        self.caching_levels.update(common_matcher.caching_levels)


class CommonParser(object):
    """
    This class is used as a base class for objects that store and instantiate
    common parts of the hierarchical SimpleMatcher structure. The object can be
    shared for many MainHierarchicalParsers.
    """
    def __init__(self, parser_context):

        # Repeating regex definitions
        self.parser_context = parser_context
        self.backend = parser_context.caching_backend
        self.file_service = parser_context.file_service
        self.cache_service = parser_context.cache_service
        self.caching_levels = {}
        self.regexs = RegexService()

    def getOnCloseTriggers(self):
        """
        Returns:
            A dictionary containing a section name as a key, and a list of
            trigger functions associated with closing that section.
        """
        onClose = {}
        for attr, callback in extractOnCloseTriggers(self).items():
            onClose[attr] = [callback]
        return onClose

    def getOnOpenTriggers(self):
        """
        Returns:
            A dictionary containing a section name as a key, and a list of
            trigger functions associated with opening that section.
        """
        onOpen = {}
        for attr, callback in extractOnOpenTriggers(self).items():
            onOpen[attr] = [callback]
        return onOpen


class ParserContext(object):
    """A container class for storing and moving information about the parsing
    environment. A single ParserContext object is initialized by the Parser
    class, or it's subclass.
    """
    def __init__(self):
        self.main_file = None
        self.version_id = None
        self.metainfo_to_keep = None
        self.super_backend = None
        self.caching_backend = None
        self.default_units = None
        self.metainfo_units = None
        self.file_service = None
        self.metainfo_env = None
        self.parser_info = None
        self.cache_service = None


class CacheObject(object):
    """Wraps a value stored inside a CacheService.
    """
    def __init__(self, name, default_value=None, single=True, update=True):
        self.name = name
        self.update = update
        self.default_value = copy.copy(default_value)
        self._value = default_value
        self._single = single
        self._pushed = False
        self._updated = default_value is not None

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value):
        self._value = value

    def clear(self):
        self._updated = True
        self._value = self.default_value


class CacheService(object):
    """A class that can be used to store intermediate results in the parsing
    process. This is quite similar to the caching system that is used by the
    ActiveBackend object but this version also allows caching values that have
    no metadata associated with them. Also you can control how many times the
    data can be read and written from the cache.
    """
    def __init__(self, parser_context=None):
        self._cache = {}
        self.parser_context = parser_context

    def __getitem__(self, name):
        """Get the value identified by name. If the cachemode does not support
        getting the value, an exception is raised.

        returns:

        raises:
        """
        cache_object = self.get_cache_object(name)
        return cache_object.value

    def get_cache_object(self, name):

        cache_object = self._cache.get(name)
        if cache_object is None:
            logger.warning("The cache object for '{}' could not be found".format(name))
            return None
        return cache_object

    def __setitem__(self, name, value):
        """Used to set the value for an item. The CacheObject corresponding to
        the name has to be first created by using the function
        add_cache_object().
        """
        cache_object = self._cache[name]
        cache_object.value = value
        cache_object._updated = True

    def clear(self):
        """Frees all object from the cache.
        """
        for cache_object in self._cache.values():
            cache_object.clear()

    def add(self, name, value=None, single=True, update=True):
        """Used to add a cache object. Two cache objects with the same name are
        not allowed.

        Args:
            single (bool): If the value is allowed to be pushed only once.
            update (bool): If the value should be update before pushing.
        """
        if name in self._cache:
            raise LookupError("There already exists a cached value with the name '{}'. All keys in the CacheService should be unique.".format(name))
        cache_object = CacheObject(name, value, single, update)
        self._cache[name] = cache_object

    def check_push_allowed(self, cache_object):
        if cache_object._single:
            if cache_object._pushed:
                raise LookupError("The CacheService value '{}' has already been output to the backend. The CacheOutputMode does not allow this.".format(cache_object.name))
                return False

        if cache_object.update:
            if not cache_object._updated:
                raise LookupError("Could not push value '{}' to backend because it was not updated since last push.".format(cache_object.name))
                return False

        return True

    def addValue(self, metaname, name=None):
        """Pushes the value stored in the cache with the given name to
        the backend.

        If the name cannot be found in the dictionary, or the value associated
        with the name is None, nothing is pushed. This allows one to use the
        cache with a simple syntax without having to worry about the actual
        availability of the data.

        If the metaname property is not defined the value is pushed with the
        name that was used as key.
        """
        if name is None:
            name = metaname

        cache_object = self.get_cache_object(name)
        if cache_object is None or cache_object._value is None:
            logger.warning("The value for metaname '{}' was not set in the CacheService, and could not be pushed".format(name))
            return

        if self.check_push_allowed(cache_object):
            self.parser_context.caching_backend.addValue(metaname, cache_object.value)
            cache_object._pushed = True
            cache_object._updated = False

    def addRealValue(self, metaname, name=None, unit=None):
        """Pushes the real value stored in the cache with the given name to
        the backend. A unit conversion will be made if units are specified.

        If the name cannot be found in the dictionary, or the value associated
        with the name is None, nothing is pushed. This allows one to use the
        cache with a simple syntax without having to worry about the actual
        availability of the data.

        If the metaname property is not defined the value is pushed with the
        name that was used as key.
        """
        if name is None:
            name = metaname

        cache_object = self.get_cache_object(name)
        if cache_object is None or cache_object._value is None:
            logger.warning("The value for metaname '{}' was not set in the CacheService, and could not be pushed".format(name))
            return

        if self.check_push_allowed(cache_object):

            if cache_object.value is not None:
                self.parser_context.caching_backend.addRealValue(metaname, cache_object.value, unit=unit)
                cache_object._pushed = True
                cache_object._updated = False
            else:
                logger.warning("The value for metaname '{}' was not set in the CacheService, and could not be pushed".format(name))

    def addArrayValues(self, metaname, name=None, unit=None):
        """Pushes the array value stored in the cache with the given name to
        the backend. A unit conversion will be made if units are specified.

        If the name cannot be found in the dictionary, or the value associated
        with the name is None, nothing is pushed. This allows one to use the
        cache with a simple syntax without having to worry about the actual
        availability of the data.

        If the metaname property is not defined the value is pushed with the
        name that was used as key.
        """
        if name is None:
            name = metaname

        cache_object = self.get_cache_object(name)
        if cache_object is None or cache_object._value is None:
            logger.warning("The value for metaname '{}' was not set in the CacheService, and could not be pushed".format(name))
            return

        if self.check_push_allowed(cache_object):

            if cache_object.value is not None:
                self.parser_context.caching_backend.addArrayValues(metaname, np.array(cache_object.value), unit=unit)
                cache_object._pushed = True
                cache_object._updated = False
            else:
                logger.warning("The value for metaname '{}' was not set in the CacheService, and could not be pushed".format(name))
