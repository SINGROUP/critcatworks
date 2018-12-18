# Copyright 2015-2018 Lauri Himanen, Fawzi Mohamed, Ankit Kariryaa
# 
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
This is the access point to the parser for the scala layer in the
NOMAD project.
"""
from __future__ import absolute_import
import sys
import setup_paths
from nomadcore.parser_backend import JsonParseEventsWriterBackend
from cp2kparser import CP2KParser


if __name__ == "__main__":

    # Initialise the parser with the main filename and a JSON backend
    main_file = sys.argv[1]
    parser = CP2KParser(backend=JsonParseEventsWriterBackend)
    parser.parse(main_file)
