echo "Be sure you are in the critcatworks directory when starting this."
echo "Do you have a new python3 virtual environment switched on?"
echo ""
echo "Automatically installing NOMAD parser, pycp2k, and critcatworks"

echo "Installing NOMAD parser"
echo "Installing NOMAD parser common"
cd nomad_parser/python-common
pip install -r requirements.txt
pip install -e .

echo "Installing NOMAD parser cp2k"
cd ../parser-cp2k
pip install -e .

echo "pycp2k - input parser installation"
echo "Installation under critcatworks/nomad_parser"

cd ..
git clone https://github.com/SINGROUP/pycp2k.git
cd pycp2k
python3 setup.py install
2
../cp2k4-1.xml
3
1

echo "Installing critcatworks"

cd ../../critcatworks

python3 setup.py install

echo "Done installing"
