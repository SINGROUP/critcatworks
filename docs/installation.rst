Installation
============

installation instructions will be improved soon


0. virtual environment recommended
virtualenv --python=python3 critcatenv
source critcatenv/bin/activate

1. Download this repository

2. Go inside this repository
cd critcatworks

3.  slightly modified NOMAD CP2K parser Installation

```sh
pip install numpy cython

cd nomad_parser/python-common
pip install -r requirements.txt
pip install -e .
```

```sh
cd ../parser-cp2k
pip install -e .

4. critcatworks installation

cd ../../

python3 setup.py install
```

5. Deactivate virtual environment

deactivate

--

steps 3 and 4 are automatically executed by 
```sh
./easily_install_all.sh
```
