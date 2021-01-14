#!/bin/sh

## Estimate changes in namelists between two releases
#####################################################

## Prerequisites: f90nml to parse Fortran namelists (http://f90nml.readthedocs.org)
## Launch under ./doc, edit the path to previous namelists accordingly

for nlst in ../../releases/release-3.6/DOC/Namelists/*; do

    if [ -f namelists/$( basename $nlst ) ]; then
        printf "%s: " $( basename $nlst )
        f90nml namelists/$( basename $nlst ) | awk '/=/ { print $0 }' | sort > nlst_new.txt
        f90nml                       $nlst   | awk '/=/ { print $0 }' | sort > nlst_old.txt

        if [[ $( diff -q nlst_*.txt ) ]]; then
            diff -y --suppress-common-lines nlst_*.txt | wc -l
        else
            echo 0
        fi

    fi

done | sort -k2rn

\rm nlst_*.txt
