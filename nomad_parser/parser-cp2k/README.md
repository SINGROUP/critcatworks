This is the main repository of the [NOMAD](https://www.nomad-coe.eu/) parser for
[CP2K](https://www.cp2k.org/).

# Example
```python
    from cp2kparser import CP2KParser
    import matplotlib.pyplot as mpl

    # 1. Initialize a parser with a set of default units.
    default_units = ["eV"]
    parser = CP2KParser(default_units=default_units)

    # 2. Parse a file
    path = "path/to/main.file"
    results = parser.parse(path)

    # 3. Query the results with using the id's created specifically for NOMAD.
    scf_energies = results["energy_total_scf_iteration"]
    mpl.plot(scf_energies)
    mpl.show()
```

# Installation
The code is python 2 and python 3 compatible. First download and install
the nomadcore package:

```sh
git clone https://gitlab.mpcdf.mpg.de/nomad-lab/python-common.git
cd python-common
pip install -r requirements.txt
pip install -e .
```

Then download the metainfo definitions to the same folder where the
'python-common' repository was cloned:

```sh
git clone https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info.git
```

Finally download and install the parser:

```sh
git clone https://gitlab.mpcdf.mpg.de/nomad-lab/parser-cp2k.git
cd parser-cp2k
pip install -e .
```

# Notes
The parser is based on CP2K 2.6.2.

The CP2K input setting
[PRINT_LEVEL](https://manual.cp2k.org/trunk/CP2K_INPUT/GLOBAL.html#PRINT_LEVEL)
controls the amount of details that are outputted during the calculation. The
higher this setting is, the more can be parsed from the upload.

The parser will try to find the paths to all the input and output files, but if
they are located very deep inside some folder structure or outside the folder
where the output file is, the parser will not be able to locate them. For this
reason it is recommended to keep the upload structure as flat as possible.

Here is a list of features/fixes that would make the parsing of CP2K results
easier:
 - The pdb trajectory output doesn't seem to conform to the actual standard as
   the different configurations are separated by the END keyword which is
   supposed to be written only once in the file. The [format
   specification](http://www.wwpdb.org/documentation/file-format) states that
   different configurations should start with MODEL and end with ENDMDL tags.
 - The output file should contain the paths/filenames of different input and
   output files that are accessed during the program run. This data is already
   available for some files (input file, most files produced by MD), but many
   are not mentioned.
