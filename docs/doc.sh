sphinx-apidoc -f -o ./src/doc ../critcatworks
cp -r _build/html/* .
make html
