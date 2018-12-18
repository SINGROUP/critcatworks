import ase.io
import ase.io.formats
import mdtraj as md
import mdtraj.formats
import numpy as np
import logging
logger = logging.getLogger("nomad")


#===============================================================================
def iread(filename, file_format=None):
    """Generator function that is used to read an atomic configuration file (MD
    trajectory, geometry optimization, static snapshot) from a file one frame
    at a time. Only the xyz positions are returned from the file, and no unit
    conversion is done, so you have to be careful with units.

    By using a generator pattern we can avoid loading the entire trajectory
    file into memory. This function will instead load a chunk of the file into
    memory (with MDTraj you can decide the chunk size, with ASE it seems to
    always be one frame), and serve individual files from that chunk. Once the
    frames in one chunk are iterated, the chunk will be garbage collected and
    memory is freed.

    Args:
        filename: String for the file path.
        file_format: String for the file format. If not given the format is
            automatically detected from the extension.

    Yields:
        numpy array containing the atomic positions in one frame.

    """
    # If file format is not explicitly stated, determine the format from the
    # filename
    if file_format is None:
        file_format = filename.split(".")[-1]

    # Try to open the file with MDTraj first. With a brief inspection it seems
    # that MDTraj is better performance wise, because it can iteratively load a
    # "chunk" of frames, and still serve the individual frames one by one. ASE
    # on the other hand will iteratively read frames one by one (unnecessary
    # IO).
    mdtraj_chunk = 100  # How many frames MDTraj will load at once
    mdtraj_failed = False

    # Must use the low level MDTraj API to open files without topology.
    class_format_map = {
            "dcd": mdtraj.formats.DCDTrajectoryFile,
            "xyz": mdtraj.formats.XYZTrajectoryFile,
            "pdb": mdtraj.formats.PDBTrajectoryFile,
    }
    traj_class = class_format_map.get(file_format)
    if traj_class is not None:
        try:
            with traj_class(filename, mode="r") as f:
                empty = False
                while not empty:
                    data = f.read(mdtraj_chunk)
                    if isinstance(data, tuple):
                        positions = data[0]
                    else:
                        positions = data
                    if len(positions) == 0:
                        empty = True
                    else:
                        for pos in positions:
                            yield pos
        except IOError:
            logger.warning("MDTraj could not read the file '{}' with format '{}'. The contents might be malformed or wrong format used.".format(filename, file_format))
            return
    else:
        mdtraj_failed = True

    # If MDTraj didn't support the format, try ASE instead
    if mdtraj_failed:
        try:
            io = ase.io.formats.get_ioformat(file_format)
        except ValueError:
            logger.error("MDTraj could not read the file '{}' with format '{}'. If MDTraj is supposed to read this format, the contents might be malformed.".format(filename, file_format))
            return
        else:
            # Return the positions in a numpy array instead of an ASE Atoms object
            generator = ase.io.iread(filename, format=file_format)
            for atoms in generator:
                pos = atoms.positions
                yield pos
