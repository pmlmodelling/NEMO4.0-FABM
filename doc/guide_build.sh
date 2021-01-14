#!/bin/sh

## Avoid the use of shell builtin echo (for -e option)
alias echo='/bin/echo -e'

## Check dependancies
##-------------------

## Sphinx, BibTeX extension and "Read The Docs" theme
if [ -n "$( ./tools/check_pkg.py sphinx sphinxcontrib.bibtex sphinx_rtd_theme )" ]; then
    echo 'One of the Python dependencies is missing => QUIT'
    exit 2
fi

cd rst

echo "\t¤ Clean previous build"
make clean
echo

echo "\t¤ Generation of the guide"
make html > /dev/null
echo

echo "\t¤ End of building run"
echo "Open ./rst/build/html/NEMO_guide.html"
cd - > /dev/null

exit 0
