from __future__ import with_statement
from __future__ import division
from __future__ import absolute_import
import os
import math
import json
import string
import h5py
import numpy as np
from abc import ABCMeta, abstractmethod
from io import open
import re

import logging


LOGGER = logging.getLogger(__name__)


class ArchiveSection(object):
    __metaclass__ = ABCMeta
    """Defines a storage independent, dictionary like interface to a section
    inside an archive file with the possibility to do recursive searches and
    indexing.

    To make a storage specific implementation just subclass this and
    define the required methods.
    """
    def __init__(self, data, path, archive):
        """
        Args:
            data: The original data that is wrapped by this ArchiveSection.
            path: The full path to this section.
            archive: The Archive object to which this section belongs to.
        """
        self._data = data
        self._path = path
        self._archive = archive

    def is_filtered(self, name):
        """Used to filter out unnecessary information when recursing the
        information. This unnecessary information includes e.g. the gIndex,
        references, etc.
        """
        filtered = set([
            "gIndex",
            "name",
            "references",
            "type",
        ])
        if name not in filtered:
            return name
        return None

    def get_by_path(self, path):
        """Used to query the ArchiveSection recursively with a simple syntax
        that also allows indexing.

        Examples:
            calculation = archive.get_by_path(
                "{repository_hash}/"
                "{calculation_hash}/"
            )
            cell = calculation.get_by_path(
                "section_run:0/"
                "section_system:0/"
                "simulation_cell"
            )
            section_method = calculation.get_by_path(
                "section_run:0/"
                "section_method:0"
            )
            multiple_sections = calculation.get_by_path(
                "section_run:0/"
                "section_method"
            )

        Args:
            path (string): A string indicating the location of the wanted
                information within the metainfo hierarchy. See the examples.

        Returns:
            One of the following: ArchiveSection, a list of ArchiveSections or a concrete value
            corresponding to the given query,
        """
        parts, _, _ = self.get_path_parts(path)
        current_data = self
        current_path = self._path
        n_parts = len(parts)
        for i_part, part in enumerate(parts):

            # Check that the section has not been deleted
            current_path = "{}/{}".format(current_path, part)
            deleted = self.check_deletions(current_path)
            if deleted:
                raise KeyError("Value for '{}' could not be found".format(current_path))

            current_data = current_data.get_child(part)
            if i_part == n_parts - 1:

                filtered = self.is_filtered(part)
                if filtered is None or filtered != part:
                    raise KeyError("Value for '{}' could not be found".format(current_path))
                return current_data

    def get(self, key, default=None):
        """Try to get a value with the given key, if none found, return the
        given default value.

        Args:
            key (string): The key to look for
            default: The value to return if the key is not found.

        Returns:
            The value for the key if found, otherwise returns the given default
            value.
        """
        # Check that the value has not been deleted
        full_path = "{}/{}".format(self._path, key)
        deleted = self.check_deletions(full_path)
        if deleted:
            return default

        # Check cache
        try:
            value = self.get_from_cache(key)
        except KeyError:
            try:
                value = self.get_by_path(key)
            except KeyError:
                value = default
        return value

    def get_scalar(self, key):
        """Get a single value from within this section. If the requested value
        is not scalar, this will raise a KeyError.

        Args:
            key (string): The key to look for

        Returns:
            A scalar value correspondin to the given key.

        Raises:
            KeyError if the found value is not sized like a scalar.
        """
        value = self[key]
        shape = value.shape
        if shape != (1,):
            raise KeyError(
                "The value for '{}/{}' is not shaped like a scalar value. "
                "Instead it has a shape of '{}'."
                .format(self._path, key, shape)
            )
        return value[0]

    def __setitem__(self, path, value):
        """Used to override or set invalid or missing values in the original
        archive.

        This method will store these values into a temporary location
        dictionary that is separate from the original data.
        """
        if self._archive.use_write_cache:

            # If this key had been deleted in the local cache, remove that mark.
            full_path = "{}/{}".format(self._path, path)
            self._archive._deletions.pop(full_path, None)

            # If the value is not already a numpy array, make it so
            if not isinstance(value, np.ndarray):
                if isinstance(value, (list, tuple)):
                    value = np.array(value)
                else:
                    value = np.array((value,))
            self._archive._overrides[full_path] = value
        else:
            raise ValueError(
                "Writing to the source file is currently disabled. If you want "
                "to write to a local cache, set the 'use_write_cache' attribute"
                " of this Archive to True."
            )

    def __delitem__(self, key):
        if self._archive.use_write_cache:
            if key in self:
                full_path = "{}/{}".format(self._path, key)
                self._archive._deletions[full_path] = True
            else:
                raise KeyError(
                    "A value for the given key '{}' has not been set, so it "
                    "could not be deleted."
                    .format(key)
                )
        else:
            raise ValueError(
                "Deleting from the source file is currently disabled. If you "
                "want to enable a local write cache, set the 'use_write_cache' "
                "attribute of this Archive to True."
            )

    def __getitem__(self, key):
        """Used to get a direct child of this section by name.

        Args:
            key (string): The name of the child entry.

        Returns:
            ArchiveSection or a concrete value.
        """
        # Check that the value has not been deleted
        full_path = "{}/{}".format(self._path, key)
        deleted = self.check_deletions(full_path)

        if deleted:
            raise KeyError("Value for '{}' has not been set.".format(full_path))
        try:
            value = self.get_from_cache(key)
        except KeyError:
            value = self.get_by_path(key)
        return value

    @abstractmethod
    def get_child(self, key):
        """This method is used to query the direct children of the section.

        This function will be called internally by the get_by_path(), get() and
        __getitem__() functions.
        """

    @abstractmethod
    def __len__(self):
        """The number of direct children.

        Returns:
            int: The number of direct children.
        """

    @abstractmethod
    def __contains__(self, key):
        """Test for whether this section contains a value for the given key.
        The key can also use the more complex syntax specified in get_by_path.

        Args:
            key (string): The key to check.

        Returns:
            boolean: Whether the given key has a value associated within this
            section.
        """

    @abstractmethod
    def items(self):
        """Used to return the keys and values of the direct children of this
        section.

        Returns:
            An iterable containing the keys and values of the direct values of
            this section.
        """

    @abstractmethod
    def keys(self):
        """Used to return the keys of the direct children of this section.

        Returns:
            An iterable containing the keys of the direct values of this
            section.
        """

    @abstractmethod
    def values(self):
        """Used to return the direct children of this section.

        Returns:
            An iterable containing the direct values of this section.
        """

    def get_from_cache(self, key):
        """Used to get a value from an internal cache that is defined for the
        Arcvive object to which this section belongs to.

        Args:
            key (string): The key to look for in the cache. The given key is
                relative to the path of this section.
        Returns:
            The value found from the cache, is it is set.
        """
        key = "{}/{}".format(self._path, key)
        value = self._archive._overrides[key]
        return value

    def check_deletions(self, full_path):
        """Used to get a value from an internal cache that is defined for the
        Arcvive object to which this section belongs to.

        Args:
            key (string): The key to look for in the cache. The given key is
                relative to the path of this section.
        Returns:
            The value found from the cache, is it is set.
        """
        value = self._archive._deletions.get(full_path)
        if value is True:
            return True
        return False

    def get_path_parts(self, path):
        """Used to separate the given path into sections, and the section to
        name and an index if specified.

        Args:
            path (string):

        Returns:
            tuple: List of section names containing also the index, a list of
                section names without the index, a list of section indices where
                the index is None if no index was specified.

        Raises:
            ValueError if the queried string has an invalid syntax.
        """
        # Remove trailing and leading slash
        if path.endswith("/"):
            path = path[:-1]
        if path.startswith("/"):
            path = path[1:]

        # Go through the hierarchy
        parts = path.split("/")
        indices = []
        names = []
        n_nones = 0
        for part in parts:
            splitted = part.split(":")
            if len(splitted) == 2:
                names.append(splitted[0])
                indices.append(int(splitted[1]))
            elif len(splitted) == 1:
                names.append(splitted[0])
                indices.append(None)
                n_nones += 1
            else:
                raise ValueError(
                    "The given path '{}' is invalid, as there are too many "
                    "colons within a section name."
                    .format(path)
                )
        if n_nones > 3:
            raise ValueError(
                "The given path contains a section without an index in the "
                "middle of the query. You will have to specify the index of "
                "all sections except the last one."
            )
        return parts, names, indices


class Archive(object):
    __metaclass__ = ABCMeta
    """Defines a storage independent interface to an archive file. To make a
    storage specific implementation just subclass this and define the required
    methods.

    Currently only reading the original Archive file is supported. You can also
    set values, but these values will not be written to the the original file
    but will only be valid during the lifetime of this object.

    Attributes:
        filepath (string): Path to the archive file
        repositories (dict): Contains all the repositories within this Archive.
        calculations (dict): Contains all the calculations within this Archive.
        use_write_cache (boolean): Whether to enable writing to a cache that
            will not persist to the original source file.
        _overrides (dictionary): Contains the values that are set during the
            lifetime of this object. These values will not persists on the
            original file.
        _deletions (dictionary): Contains the keys that are deleted during the
            lifetime of this object. These deletions will not be made to the
            original file.
    """
    def __init__(self, filepath, use_write_cache=False):
        """
        Args:
            filepath (string): Filepath to an archive file.
            use_write_cache (boolean): Whether to enable writing to a cache
                that will not persist to the original source file.
        """
        self.filepath = filepath
        self.repositories = {}
        self.calculations = {}
        self.use_write_cache = use_write_cache
        self._overrides = {}
        self._deletions = {}

    @staticmethod
    def factory(archive_path, use_write_cache=False):
        """A factory method for creating Archive objects based on the file
        type.

        Args:
            archive_path (string): The path to the archive file.

        Returns:
            An Archive subclass that is specific to the filetype of the given
            file.
        """
        extension = archive_path.rsplit(".", 1)[-1]
        if extension == "json":
            return ArchiveJSON(archive_path, use_write_cache)
        elif extension == "h5":
            return ArchiveHDF5(archive_path, use_write_cache)
        else:
            raise ValueError(
                "Unknown archive filetype with extension '{}'."
                .format(extension)
            )

    @abstractmethod
    def setup(self, root):
        """Used to setup the dictionaries that contain the repositories and
        calculations.
        """


class ArchiveHDF5(Archive):
    """Provides an access to a HDF5 archive files.

    Also keeps a local cache for the index data. The same index data is often
    queried multiple times, so caching it improves the access speed. The index
    data is also so small that keeping a local cache of it doesn't consume much
    memory.

    Attributes:
        index_cache (dict): A cache containing the index data for groups and
            datasets.
    """
    def __init__(self, filepath, use_write_cache=False):
        super(ArchiveHDF5, self).__init__(filepath, use_write_cache)
        h5_root = h5py.File(filepath, "r")
        self.index_cache = {}
        self.setup(h5_root)

    def setup(self, root):
        """Used to setup the dictionaries that contain the repositories and
        calculations.
        """
        for repo_name, repo in root.items():
            self.repositories[repo_name] = {}
            for calc_name, calc in repo.items():
                calc_section = self.calculations[calc_name] = ArchiveSectionHDF5(
                    calc,
                    "{}/{}".format(repo_name, calc_name),
                    self,
                    [[0]],
                    0
                )
                self.calculations[calc_name] = calc_section
                self.repositories[repo_name][calc_name] = calc_section


class ArchiveSectionHDF5(ArchiveSection):
    """A section from a HDF5 file. Roughly corresponds to the concept of a
    Group in HDF5.
    """
    BASE64DIGITS = string.ascii_uppercase + string.ascii_lowercase + string.digits + "+" + "/"

    def __init__(self, data, path, archive, index_datas, local_index):
        super(ArchiveSectionHDF5, self).__init__(data, path, archive)
        _, names, indices = self.get_path_parts(path)
        # Here we drop out the indices of the repository and calculation
        # section, as they are "None"
        self._indices = indices[2:]
        self._names = names
        self._index_datas = index_datas
        self._local_index = local_index

    def __len__(self):
        return len(self.keys())

    def is_filtered(self, key):
        # return key
        if key.endswith("-index"):
            return None
        key_without_size = key.rsplit(".", 1)[0]
        if key_without_size.endswith("-v"):
            return key_without_size[:-2]
        else:
            return key

    def __contains__(self, key):
        try:
            self[key]
        except KeyError:
            return False
        return True

    def items(self):
        local_index = 0
        for key, value in self._data.items():
            key_filtered = self.is_filtered(key)
            if key_filtered is None:
                continue
            if isinstance(value, h5py.Group):
                index_datas = self._index_datas[:]
                index_data = self._data.get("{}-index".format(key))
                if index_data is None:
                    index_datas = []
                else:
                    index_datas.append(index_data)
                yield (key_filtered, ArchiveSectionHDF5(
                    value,
                    "{}/{}".format(self._path, key),
                    self._archive,
                    index_datas,
                    local_index)
                )
            else:
                yield (key_filtered, value)
            local_index += 1

    def keys(self):
        for key in self._data.keys():
            key_filtered = self.is_filtered(key)
            if key_filtered is not None:
                yield key_filtered

    def values(self):
        local_index = 0
        for key, value in self._data.items():
            key_filtered = self.is_filtered(key)
            if key_filtered is None:
                continue
            if isinstance(value, h5py.Group):
                index_datas = self._index_datas[:]
                index_data = self._data.get("{}-index".format(key))
                if index_data is None:
                    index_datas = []
                else:
                    index_datas.append(index_data)
                yield ArchiveSectionHDF5(
                    value,
                    "{}/{}".format(self._path, key),
                    self._archive,
                    index_datas,
                    local_index
                )
            else:
                yield value
            local_index += 1

    @property
    def mainfile_uri(self):
        """Used to the the mainfile uri from a section corresponding to a
        calculation.
        """
        numpy_uri = self._data.attrs["mainFileUri"]
        str_uri = numpy_uri[0].decode("utf-8")
        return str_uri

    @property
    def parser_info(self):
        """Used to the parser info from a section corresponding to a
        calculation.
        """
        numpy_uri = self._data.attrs["parserInfo"]
        parser_info_str = numpy_uri[0].decode("utf-8")
        parser_info_dict = json.loads(parser_info_str)
        return parser_info_dict

    def check_indices(self):
        """This function is used to check for checking the
        """

    def get_child(self, path):
        """Used to get a direct child of this section within the HDF5 file.
        """
        # Full name for the object that is being searched.
        child_path = "{}/{}".format(self._path, path)

        # Create an index variable
        child_index = None
        if len(self._indices) == 0:
            indices = [0]
        else:
            indices = self._indices[:]
        parent_index = indices[-1]

        # Create a copy of the index data for the child section that is
        # retrieved.
        index_datas = self._index_datas[:]

        # Separate the section number and name
        splitted = path.rsplit(":", 1)
        if len(splitted) == 2:
            name = splitted[0]
            child_index = int(splitted[1])
            indices.append(child_index)
        elif len(splitted) == 1:
            name = path
            child_index = None

        # Check if the path points to a section or a value. Sections will be
        # found simply with the metaname, but the concrete values have a -v and
        # possibly the shape appended to them
        section = True
        try:
            data = self._data[name]
        except KeyError:
            section = False

        # Get the new index data for this child. First try to get the data from
        # a local cache, if not present, then look up from the HDF5 file.
        if section:
            index_path = "{}/{}-index".format(name, name)
        else:
            index_path = "{}-index".format(name)
        global_index_path = "/".join(self._names) + "/" + index_path
        index_data = self._archive.index_cache.get(global_index_path)
        if index_data is None:
            index_data = self._data.get(index_path)
            if index_data is None:
                index_data = np.array([[0]])
            else:
                index_data = index_data.value
            self._archive.index_cache[global_index_path] = index_data
        index_datas.append(index_data)

        # If a section is requested, either provide the section with the
        # provided index, or then a list of section if index was not specified.
        if section:
            # Get the list of available parent indices
            parent_indices = index_data[:, 0]

            # See which section in the index list have the correct parent,
            # and get the correct child from there
            candidates = np.where(parent_indices == parent_index)[0]
            n_candidates = len(candidates)
            if n_candidates == 0:
                raise KeyError(
                    "Could not find value at path '{}'."
                    .format(child_path)
                )
            if child_index is not None:
                try:
                    child_index = candidates[child_index]
                except IndexError:
                    raise KeyError(
                        "Could not find value at path '{}'."
                        .format(child_path)
                    )

                section = self._data[name]
                data = ArchiveSectionHDF5(
                    section,
                    child_path,
                    self._archive,
                    index_datas,
                    child_index
                )
            else:
                group_list = []
                section = self._data[name]
                for candidate in candidates:

                    # Only add the sections that have not been deleted
                    full_path = "{}/{}:{}".format(self._path, name, candidate)
                    deleted = self.check_deletions(full_path)
                    if not deleted:
                        group = ArchiveSectionHDF5(
                            section,
                            full_path,
                            self._archive,
                            index_datas,
                            candidate
                        )
                        group_list.append(group)
                data = group_list

        # If a dataset is requested, check that the dataset exists for this
        # path and then provide it as a numpy array.
        else:
            # Get the local indices of the parent section and search for them
            # in the first column of the index table. There should only be one
            # entry. The index of the data associated with this entry is
            # indicated by the second column. The dataset name can be decided
            # from the data in the remaining columns, which indicate the
            # dimensions of the data.
            test_index = np.where(index_data[:, 0] == self._local_index)[0]
            if test_index.size == 0:
                raise KeyError(
                    "Could not find value at path '{}'."
                    .format(child_path)
                )

            # This error is currently disabled, because it seems that the
            # metainfo system supports repeating scalar values for one section.
            # if test_index.size > 1:
                # raise ValueError(
                    # "The HDF file contains more than one dataset for the "
                    # "path '{}'. "
                    # .format(child_path)
                # )

            # If the value can have multiple shapes, the values are split into
            # different tables. For each table there is a local index in the
            # second column of the index table that we must use.
            data = []
            for row_i in test_index:
                index_row = index_data[row_i]
                if index_row.shape != (1,):
                    data_index = index_row[1]
                else:
                    data_index = row_i

                # The data name may depend on the shape, and if so, the
                # shape is appended to the name as base64 fields
                data_path = name + "-v"
                index_shape = index_data.shape
                if index_shape[1] > 2:
                    data_path = name + "-v"
                    for dim in index_data[data_index][2:]:
                        base64dim = self.base64convert(dim)
                        data_path += ".{}".format(base64dim)

                i_data = self._data[data_path][data_index]

                # Convert bytestrings to regular strings
                if i_data.dtype == np.object:
                    i_data = np.array([x.decode("utf-8") for x in i_data])

                # Gather scalar values to a 1D list
                if i_data.shape == (1,):
                    data.append(i_data[0])
                else:
                    data.append(i_data)

            # If one object returned, remove the outermost list
            if len(test_index) == 1:
                if data[0].shape == ():
                    data = np.array([data[0]])
                else:
                    data = data[0]
            else:
                data = np.array(data)

        return data

    def base64convert(self, x):
        """Used to convert the given base 10 number to base 64.
        """
        base = 64
        if x == 0:
            return ArchiveSectionHDF5.BASE64DIGITS[0]
        digits = []

        while x:
            digits.append(ArchiveSectionHDF5.BASE64DIGITS[x % base])
            x = math.floor(x/base)

        return "".join(digits)


class ArchiveJSON(Archive):
    """Provides a dicionary like access to JSON archive files with the
    possibility to do recursive searches and indexing.

    This implementation will load the entire JSON file into memory, which might
    become a problem with big files and parallel execution on the same machine.
    """
    def __init__(self, filepath, use_write_cache=False):
        super(ArchiveJSON, self).__init__(filepath, use_write_cache)
        with open(filepath, "r") as fin:

            json_root = json.load(fin)

            # Get the calculation name from filename
            calculation_name = os.path.basename(filepath).split(".", 1)[0]

            # Get the repository name from mainFileUri
            mainfile_uri = json_root["mainFileUri"]
            repository_name = mainfile_uri.split("://", 1)[1]
            repository_name = repository_name.split("/", 1)[0]

            root_section = {
                repository_name: {
                    calculation_name: json_root
                }
            }
            self.setup(root_section)

    def setup(self, root):
        for repo_name, repo in root.items():
            self.repositories[repo_name] = {}
            for calc_name, calc in repo.items():
                calc_section = self.calculations[calc_name] = ArchiveSectionJSON(
                    calc,
                    "{}/{}".format(repo_name, calc_name),
                    self
                )
                self.calculations[calc_name] = calc_section
                self.repositories[repo_name][calc_name] = calc_section


class ArchiveSectionJSON(ArchiveSection):
    """Represents a section inside a JSON-file.
    """
    def __len__(self):
        return len(self._data)

    def __contains__(self, key):
        try:
            self[key]
        except KeyError:
            return False
        return True

    def items(self):
        for key, value in self._data.items():
            if self.is_filtered(key) is None:
                continue
            if isinstance(value, dict):
                yield (key, ArchiveSectionJSON(value, "{}/{}".format(self._path, key), self._archive))
            else:
                # Wrap scalar values inside a numpy array to be consistent with the HDF5
                # Archive.
                if not isinstance(value, list):
                    value = np.array([value])
                else:
                    value = np.array(value)
                yield (key, value)

    def keys(self):
        for key in self._data.keys():
            if self.is_filtered(key) is None:
                continue
            yield key

    def values(self):
        for key, value in self._data.items():
            if self.is_filtered(key) is None:
                continue
            if isinstance(value, dict):
                yield ArchiveSectionJSON(value, "{}/{}".format(self._path, key), self._archive)
            else:
                # Wrap scalar values inside a numpy array to be consistent with the HDF5
                # Archive.
                if not isinstance(value, list):
                    value = np.array([value])
                else:
                    value = np.array(value)
                yield value

    @property
    def mainfile_uri(self):
        """Used to the the mainfile uri from a section corresponding to a
        calculation.
        """
        return self._data["mainFileUri"]

    @property
    def parser_info(self):
        """Used to the the parser info from a section corresponding to a
        calculation.
        """
        return self._data["parserInfo"]

    def get_child(self, path):

        def get_root_sections(name):
            section_found = False
            sections = []
            i_section = 0
            while not section_found:
                try:
                    data = self._data["sections"]["{}-{}".format(name, i_section)]
                    sections.append(data)
                except KeyError:
                    section_found = True
                i_section += 1
            if len(sections) == 0:
                raise KeyError("No sections found within the path '{}'".format(path))
            return sections

        # Full name for the object that is being searched.
        child_path = "{}/{}".format(self._path, path)

        # Separate the name and index
        splitted = path.rsplit(":", 1)
        if len(splitted) == 2:
            name = splitted[0]
            index = int(splitted[1])
        else:
            name = splitted[0]
            index = None

        # Determine if section requested
        is_section = False
        if path.startswith("section"):
            is_section = True
        elif re.match(r'^x_\S+_section', path):
            # code-specific section
            is_section = True

        # If no index specified, try to get as concrete value or as a list of
        # sections
        if index is None:
            if not is_section:
                data = self._data[name]
                # Wrap scalar values inside a numpy array to be consistent with the HDF5
                # Archive.
                if not isinstance(data, list):
                    data = np.array([data])
                else:
                    data = np.array(data)
            else:
                try:
                    sections = self._data[name]
                # If this section is section_run, the children will have the
                # index in the section name directly (only specific to JSON)
                except KeyError:
                    sections = get_root_sections(name)

                data = []

                # Only add the sections that have not been deleted
                for i_section, sec in enumerate(sections):
                    full_path = "{}:{}".format(child_path, i_section)
                    deleted = self.check_deletions(full_path)
                    if not deleted:
                        data.append(ArchiveSectionJSON(sec, full_path, self._archive))

        # If index specified, try to get the specific section.
        elif is_section:
            # The section_run and it's child section are within a separate
            # "sections" dictionary, and have indices directly in the name
            try:
                data = self._data["sections"]["{}-{}".format(name, index)]
            except KeyError:
                data = self._data[name]
                try:
                    data = data[index]
                except IndexError:
                    raise KeyError("No sections found within the path '{}'".format(path))
            data = ArchiveSectionJSON(data, child_path, self._archive)
        else:
            raise ValueError(
                "Indexing within a concrete value is not supported. You "
                "will first have to get the whole data and index it then."
            )

        return data
