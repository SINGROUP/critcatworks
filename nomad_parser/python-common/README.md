# __python-common__

This repository contains the common python files that are
part of the [NOMAD Laboratory](http://nomad-lab.eu).
The official version lives at

    git@gitlab.mpcdf.mpg.de:nomad-lab/python-common.git

you can browse it at

    https://gitlab.mpcdf.mpg.de/nomad-lab/python-common

Some things rely on having the nomad-meta-info checked out at the same level.
The simplest way to have this is to check out nomad-lab-base recursively:

    git clone --recursive git@gitlab.mpcdf.mpg.de:nomad-lab/nomad-lab-base.git

then this will be in python-common within it.

# Installation
The code is python>=2.7 and python>=3.4 compatible. First download and install
this package:

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
