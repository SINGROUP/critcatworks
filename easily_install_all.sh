echo "Be sure you are in the critcatworks directory when starting this."
echo "Do you have a new python3 virtual environment switched on?"
echo ""
echo "Automatically installing NOMAD parser, and critcatworks"

pip install numpy cython

echo "Installing NOMAD parser"
echo "Installing NOMAD parser common"
cd nomad_parser/python-common
pip install -r requirements.txt
pip install -e .

echo "Installing NOMAD parser cp2k"
cd ../parser-cp2k
pip install -e .

echo "Installation under critcatworks/nomad_parser"

echo "Installing critcatworks"

cd ../..

pip install -e .

echo "Done installing"
