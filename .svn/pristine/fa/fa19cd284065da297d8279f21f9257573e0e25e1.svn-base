#!/bin/sh

## Initialisation
##---------------

## Default selection for models
if   [ "$1" = 'all' ]; then
    models='NEMO SI3 TOP'
elif [ "$1" =    '' ]; then
    models='NEMO'
else
    models=$*
fi

# Source shared functions
. tools/shr_func.sh


## Check dependancies
##-------------------

## LaTeX installation, find latexmk should be enough
[ -z $( which latexmk ) ] && { echo 'latexmk not installed => QUIT'; exit 2; }

## Pygments package for syntax highlighting of source code (namelists & snippets)
[ -n "$( ./tools/check_pkg.py pygments )" ] && { echo 'Python pygments is missing => QUIT'; exit 2; }

## Retrieve figures if not already there
#if [ ! -d latex/figures ]; then
#    printf "Downloading of shared figures and logos\n\n"
#    svn co http://forge.ipsl.jussieu.fr/nemo/svn/utils/figures latex/figures > /dev/null
#fi


## Loop on the models
##-------------------

for model in $models; do
    echo $model
#    clean $model
    build $model
    printf "\t¤ End of building run\n"
    echo
done

exit 0
