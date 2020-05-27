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
Setups the python-common library in the PYTHONPATH system variable.
"""
import sys
import os
import os.path

baseDir = os.path.dirname(os.path.abspath(__file__))
commonDir = os.path.normpath(os.path.join(baseDir, "../../../../../python-common/common/python"))
parserDir = os.path.normpath(os.path.join(baseDir, "../../parser-cp2k"))

# Using sys.path.insert(1, ...) instead of sys.path.insert(0, ...) based on
# this discusssion:
# http://stackoverflow.com/questions/10095037/why-use-sys-path-appendpath-instead-of-sys-path-insert1-path
if commonDir not in sys.path:
    sys.path.insert(1, commonDir)
if parserDir not in sys.path:
    sys.path.insert(1, parserDir)
