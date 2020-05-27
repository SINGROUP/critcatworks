import numpy as np
import logging
import re
logger = logging.getLogger(__name__)


#===============================================================================
def iread(filepath, columns=None, delimiter=r"\s+", comments=r"#", start=None, end=None, n_conf=None):
    """Used to iterate a CSV-like file. If a separator is provided the file
    is iterated one configuration at a time. Only keeps one configuration
    of the file in memory. If no separator is given, the whole file will be
    handled.

    The contents are separated into configurations whenever the separator
    regex is encountered on a line.

    Args:
        filepath: Path to the CSV like file to be processed.
        columns: List of integers indicating the columns of interest in the CSV file.
        start: A regex that is used to indicate the start of a new configuration.
        end: A regex that is used to indicate the end of a configuration.
        comments: A regex that is used identify comments in the file that are ignored.
        n_conf: Number of lines in a configuration. If you want to use multiple
            lines as a single configuration.
    """

    def split_line(line):
        """Chop off comments, strip, and split at delimiter.
        """
        if line.isspace():
            return None
        if comments:
            line = compiled_comments.split(line, maxsplit=1)[0]
        line = line.strip('\r\n ')
        if line:
            return compiled_delimiter.split(line)
        else:
            return None

    def is_end(line):
        """Check if the given line matches the separator pattern.
        Separators are used to split a file into multiple configurations.
        """
        if end:
            return compiled_end.search(line)
        return False

    def is_start(line):
        """Check if the given line matches the separator pattern.
        Separators are used to split a file into multiple configurations.
        """
        if start:
            return compiled_start.search(line)
        return False

    # Precompile the different regexs before looping
    compiled_delimiter = re.compile(delimiter)
    if comments:
        comments = (re.escape(comment) for comment in comments)
        compiled_comments = re.compile('|'.join(comments))
    if end:
        compiled_end = re.compile(end)
    if start:
        compiled_start = re.compile(start)

    # Columns as list
    if columns is not None:
        columns = list(columns)

    # Start iterating
    configuration = []
    started = False

    # If no starting and ending condition are provided, read configuration by line
    if start is None and end is None and n_conf is None:
        with open(filepath, "r") as f:
            for line in f:  # This actually reads line by line and only keeps the current line in memory
                # Ignore comments, separate by delimiter
                vals = split_line(line)
                line_forces = []
                if vals:
                    for column in columns:
                        try:
                            value = vals[column]
                        except IndexError:
                            logger.warning("The given index '{}' could not be found on the line '{}'. The given delimiter or index could be wrong.".format(column, line))
                            return
                        try:
                            value = float(value)
                        except ValueError:
                            logger.warning("Could not cast value '{}' to float. Currently only floating point values are accepted".format(value))
                            return
                        else:
                            line_forces.append(value)
                    yield np.array(line_forces)

    # If starting and ending condition are provided, after starting condition
    # is detected, add the values from lines to a new array that is returned
    # when the end condition is met
    elif start is not None and end is not None:
        with open(filepath, "r") as f:
            for line in f:  # This actually reads line by line and only keeps the current line in memory

                # If a start regex is provided, use it to detect the start of a configuration
                if is_start(line):
                    started = True
                    continue

                # If separator encountered, yield the stored configuration
                if is_end(line):
                    started = False
                    if configuration:
                        yield np.array(configuration)
                        configuration = []

                elif start is not None and started:
                    # Ignore comments, separate by delimiter
                    vals = split_line(line)
                    line_forces = []
                    if vals:
                        for column in columns:
                            try:
                                value = vals[column]
                            except IndexError:
                                logger.warning("The given index '{}' could not be found on the line '{}'. The given delimiter or index could be wrong.".format(column, line))
                                return
                            try:
                                value = float(value)
                            except ValueError:
                                logger.warning("Could not cast value '{}' to float. Currently only floating point values are accepted".format(value))
                                return
                            else:
                                line_forces.append(value)
                        configuration.append(line_forces)

            # The last configuration is yielded even if separator is not present at
            # the end of file or is not given at all
            if configuration:
                yield np.array(configuration)

    # If n_conf is defined, read multiple lines as one configuration
    elif start is None and end is None and n_conf is not None:
        with open(filepath, "r") as f:
            i_line = 0
            conf = []
            for line in f:  # This actually reads line by line and only keeps the current line in memory
                # Ignore comments, separate by delimiter
                vals = split_line(line)
                line_values = []
                if vals:
                    for column in columns:
                        try:
                            value = vals[column]
                        except IndexError:
                            logger.warning("The given index '{}' could not be found on the line '{}'. The given delimiter or index could be wrong.".format(column, line))
                            return
                        try:
                            value = float(value)
                        except ValueError:
                            logger.warning("Could not cast value '{}' to float. Currently only floating point values are accepted".format(value))
                            return
                        else:
                            line_values.append(value)
                    conf.append(line_values)
                    i_line += 1
                    if i_line == n_conf:
                        yield np.array(conf)
                        conf = []
                        i_line = 0
