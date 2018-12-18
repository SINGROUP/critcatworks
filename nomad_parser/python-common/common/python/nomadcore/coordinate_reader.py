from builtins import next
from builtins import object
import ase.io
import logging
logger = logging.getLogger(__name__)


#===============================================================================
class CoordinateReader(object):
    """Used to parse various different atomic coordinate files.

    See the dictionary 'formats' for all the supported formats and a brief
    explanation. Reading is primarily done by ASE, but in some cases own
    implementation must be used. Returns all coordinates as numpy arrays.
    """
    def __init__(self):
        self.formats = {
            "xyz":    CoordinateFormat(".xyz", "The XYZ file format.", self.ase_iread),
            "cif":    CoordinateFormat(".cif", "Crystallographic Information File", self.ase_iread),
            "pdb":    CoordinateFormat(".pdb", "Protein Data Bank", self.ase_iread),
        }

    def iread(self, file_handle, format):
        """Returns an iterator that goes through the given trajectory file one
        configuration at a time.

        Args:
            file_handle: A file object pointing to the coordinate file.
            format: String containing one of the supported formats.
        """
        if not self.check_format_support(format):
            return

        iread_function = self.formats[format].function
        return iread_function(file_handle, format)

    def n_atoms(self, file_handle, format):
        """Read the first configuration of the coordinate file to extract the
        number of atoms in it.
        """
        iterator = self.iread(file_handle, format)
        pos = next(iterator)
        return pos.shape[0]

    def check_format_support(self, format):
        """Check if the given format is supported.
        """

        if format not in self.formats:
            logger.error("The format '{}' is not supported by CoordinateReader.".format(format))
            return False
        else:
            return True

    def ase_iread(self, file_handle, format):
        """Wrapper for ASE's iread function. Returns numpy arrays instead of
        Atoms objects.
        """
        # The newest ASE version found in Github has an iread function.
        # After reading the ASE source code, it seems that the ASE iread does
        # actually read the entire file into memory and the yields the
        # configurations from it. Should be checked at some point.
        def ase_generator(iterator):
            """Used to wrap an iterator returned by ase.io.iread so that it returns
            the positions instead of the ase.Atoms object.
            """
            for value in iterator:
                yield value.get_positions()

        iterator = ase.io.read(file_handle, index=":", format=format)
        return ase_generator(iterator)

    def custom_iread(self, file_handle, format):
        """
        """
        pass


#===============================================================================
class CoordinateFormat(object):
    """Represents a coordinate format.
    """
    def __init__(self, extension, info, function):
        self.extension = extension
        self.info = info
        self.function = function
