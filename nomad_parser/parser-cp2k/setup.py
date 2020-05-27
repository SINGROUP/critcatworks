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

from setuptools import setup, find_packages


def main():
    setup(
        name="cp2kparser",
        version="0.1",
        package_data={
            'cp2kparser.versions.cp2k262': ['input_data/cp2k_input_tree.pickle'],
        },
        description="NOMAD parser implementation for CP2K.",
        author="Lauri Himanen",
        author_email="lauri.himanen@aalto.fi",
        license="GPL3",
        package_dir={'': 'parser/parser-cp2k'},
        packages=find_packages(),
        install_requires=[
            'nomadcore',
        ],
    )

if __name__ == "__main__":
    main()
