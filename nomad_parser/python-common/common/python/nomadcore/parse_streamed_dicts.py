from future import standard_library
standard_library.install_aliases()
from builtins import object
import json
import io

# the 3 possible states while reading input character per character
IN_NORMAL_TEXT = 0
IN_STRING = 1
IN_STRING_ESCAPE = 2

def readDict(inF, line0):
    "reads a dictionary from an indented json"
    status = IN_NORMAL_TEXT
    # iterator for character position in string
    i = 0
    # counts number of open {
    n_paren = 0
    # used for StringIO if dictionary is distributed over several blocks
    fullDictStr = None
    line = line0
    iline = 0
    while True:
        # iterate through characters of current block
        while i < len(line):
            c = line[i]
            i += 1
            if status == IN_STRING_ESCAPE:
                status = IN_STRING
                continue
            elif status == IN_NORMAL_TEXT:
                # skip characters that do not change state
                while (c != '{' and c!= '}' and c != '"' and i < len(line)):
                    c = line[i]
                    i += 1
                if c == '{':
                    # skip the characters before the first opening {
                    if n_paren == 0:
                        line = line[i - 1:]
                        i = 1
                    n_paren += 1
                elif c == '}':
                    n_paren -= 1
                    # found closing }, return dictionary
                    if n_paren == 0:
                        # write block upto current character to string or StringIO
                        # only use StringIO if dictionary is distributed over several blocks
                        if fullDictStr == None:
                            outS = line[:i]
                        else:
                            fullDictStr.write(line[:i])
                            outS = fullDictStr.getvalue()
                            # reset StringIO
                            fullDictStr.close()
                            fullDictStr = None
                        # reset block to remaining characters
                        line = line[i:]
                        i = 0
                        # dictionary output
                        try:
                            outD = json.loads(outS)
                        except:
                            raise Exception("Could not convert string " + repr(outS) + " with json.loads to dictionary.")
                        return outD
                elif c == '"':
                    status = IN_STRING
            elif status == IN_STRING:
                # skip characters that do not change state
                while (c != '"' and c != '\\' and i < len(line)):
                    c = line[i]
                    i += 1
                if c == '"':
                    status = IN_NORMAL_TEXT
                elif c == '\\':
                    status = IN_STRING_ESCAPE
        # if we arrive here, then the dictionary is distributed over several blocks
        # write block to StringIO but only if an opening { already found
        if n_paren > 0:
            if fullDictStr == None:
                fullDictStr = io.StringIO()
            fullDictStr.write(line)
        else:
            # failed to find dictionary in first line, stopping
            return None
        # read new block from input
        line = inF.readline()
        if not line:
            # early EOF
            return None
        iline += 1
        i = 0

def readArray(inF, line0):
    "reads an array from an indented json"
    status = IN_NORMAL_TEXT
    # iterator for character position in string
    i = 0
    # counts number of open [
    n_paren = 0
    # used for StringIO if dictionary is distributed over several blocks
    fullDictStr = None
    line = line0
    iline = 0
    while True:
        # iterate through characters of current block
        while i < len(line):
            c = line[i]
            i += 1
            if status == IN_STRING_ESCAPE:
                status = IN_STRING
                continue
            elif status == IN_NORMAL_TEXT:
                # skip characters that do not change state
                while (c != '[' and c!= ']' and c != '"' and i < len(line)):
                    c = line[i]
                    i += 1
                if c == '[':
                    # skip the characters before the first opening [
                    if n_paren == 0:
                        line = line[i - 1:]
                        i = 1
                    n_paren += 1
                elif c == ']':
                    n_paren -= 1
                    # found closing ], return dictionary
                    if n_paren == 0:
                        # write block upto current character to string or StringIO
                        # only use StringIO if dictionary is distributed over several blocks
                        if fullDictStr == None:
                            outS = line[:i]
                        else:
                            fullDictStr.write(line[:i])
                            outS = fullDictStr.getvalue()
                            # reset StringIO
                            fullDictStr.close()
                            fullDictStr = None
                        # reset block to remaining characters
                        line = line[i:]
                        i = 0
                        # dictionary output
                        try:
                            outD = json.loads(outS)
                        except:
                            raise Exception("Could not convert string " + repr(outS) + " with json.loads to dictionary.")
                        return outD
                elif c == '"':
                    status = IN_STRING
            elif status == IN_STRING:
                # skip characters that do not change state
                while (c != '"' and c != '\\' and i < len(line)):
                    c = line[i]
                    i += 1
                if c == '"':
                    status = IN_NORMAL_TEXT
                elif c == '\\':
                    status = IN_STRING_ESCAPE
        # if we arrive here, then the dictionary is distributed over several blocks
        # write block to StringIO but only if an opening [ already found
        if n_paren > 0:
            if fullDictStr == None:
                fullDictStr = io.StringIO()
            fullDictStr.write(line)
        else:
            # failed to find dictionary in first line, stopping
            return None
        # read new block from input
        line = inF.readline()
        if not line:
            # early EOF
            return None
        iline += 1
        i = 0

class ParseStreamedDicts(object):
    """allows the extraction of JSON dictionaries out of file objects
    which are then converted to python dictionaries
    therefore, strings must be passed as StringIO objects"""
    def __init__(self, inF, blockSize = 8192):
        self.inF = inF
        # read at least 2 unicode characters (2 x 4 bytes)
        if blockSize < 8:
            self.blockSize = 8
        else:
            self.blockSize = blockSize
        # read input in blocks
        self.blockRead = ""
        # set initial state
        self.status = IN_NORMAL_TEXT
        # iterator for character position in string
        self.i = 0
        # counts number of open {
        self.n_paren = 0 
        # used for StringIO if dictionary is distributed over several blocks
        self.fullDictStr = None

    def readNextDict(self):
        """reads input until a complete set of {} is found
        the so obtained string is then converted to a python dictionary with json.loads"""
        while True:
            # iterate through characters of current block
            while self.i < len(self.blockRead):
                c = self.blockRead[self.i]
                self.i += 1
                if self.status == IN_STRING_ESCAPE:
                    self.status = IN_STRING
                    continue
                elif self.status == IN_NORMAL_TEXT:
                    # skip characters that do not change state
                    while (c != '{' and c!= '}' and c != '"' and self.i < len(self.blockRead)):
                        c = self.blockRead[self.i]
                        self.i += 1
                    if c == '{':
                        # skip the characters before the first opening {
                        if self.n_paren == 0:
                            self.blockRead = self.blockRead[self.i - 1:]
                            self.i = 1
                        self.n_paren += 1
                    elif c == '}':
                        self.n_paren -= 1
                        # found closing }, return dictionary
                        if self.n_paren == 0:
                            # write block upto current character to string or StringIO
                            # only use StringIO if dictionary is distributed over several blocks
                            if self.fullDictStr == None:
                                outS = self.blockRead[:self.i] 
                            else:
                                self.fullDictStr.write(self.blockRead[:self.i])
                                outS = self.fullDictStr.getvalue()
                                # reset StringIO
                                self.fullDictStr.close()
                                self.fullDictStr = None
                            # reset block to remaining characters
                            self.blockRead = self.blockRead[self.i:]
                            self.i = 0
                            # dictionary output
                            try:
                                outD = json.loads(outS)
                            except:
                                raise Exception("Could not convert string " + repr(outS) + " with json.loads to dictionary.")
                            return outD
                    elif c == '"':
                        self.status = IN_STRING
                elif self.status == IN_STRING:
                    # skip characters that do not change state
                    while (c != '"' and c != '\\' and self.i < len(self.blockRead)):
                        c = self.blockRead[self.i]
                        self.i += 1
                    if c == '"':
                        self.status = IN_NORMAL_TEXT
                    elif c == '\\':
                        self.status = IN_STRING_ESCAPE
            # if we arrive here, then the dictionary is distributed over several blocks
            # write block to StringIO but only if an opening { already found
            if self.n_paren > 0:
                if self.fullDictStr == None:
                    self.fullDictStr = io.StringIO()
                self.fullDictStr.write(self.blockRead)
            # read new block from input
            self.blockRead = self.inF.read(self.blockSize)
            self.i = 0
            # reached end of input
            if not self.blockRead:
                return None

