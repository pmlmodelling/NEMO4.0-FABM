#! /bin/sh
#
#       usage for NEMO doc to create an update of the Namelists directory :
# 1- delete the existing directory (You can also choose to save it somewhere) 
#    rm -rf Namelists
# 2- create the updated Namelists directory from the SHARED/namelist_ref :
#    ./namelist_split.sh -i ../NEMOGCM/CONFIG/SHARED/namelist_ref -o Namelists
#
# .. _namelist_split.sh:
#
# =================
# namelist_split.sh
# =================
#
# ----------------
# split a namelist
# ----------------
#
# SYNOPSIS
# ========
#
# .. code-block:: bash
#
#    namelist_split.sh -i namelist -o dirout [-n]
#
# DESCRIPTION
# ===========
#
# Split a namelist file (NEMO convention syntax) given in parameter
# into files in a given output directory.
#
# .. option:: -i <input file (namelist)>
# .. option:: -o <output directory>
# .. option:: -n <numbered output files>
#
# Each file of the output directory is named after the *block* it contains.
#
# We assume here that the input file is written with this pattern of block:
#
# .. parsed-literal::
#
#    !---------------
#    &nameblock
#    !---------------
#    var1 = val1
#    ...
#    varn = valn
#    /
#
# If ``-n`` option is used, files in the directory output are named after
# the number of the block it contains and the name of this block. This might
# be useful to recombine blocks in the same order.
#
# EXAMPLES
# ========
#
# To split the namelist NEMO ORCA2_LIM experiment 00 if you are in NEMO
# working space :
#
# .. code-block:: bash
#
#    namelist_split.sh -i CONFIG/ORCA2_LIM/EXP00/namelist -o /tmp/EXP00_namelist_split/
#
# To split the namelist under data directory with numbered files:
#
# .. code-block:: bash
#
#    namelist_split.sh -i ./data/namelist -o /tmp/EXP00_namelist_split/ -n
#
# Check can be made by concatenation of files in output directory if ``-n``
# option have been used and comparison with original file :
#
# .. code-block:: bash
#
#    cat /tmp/EXP00_namelist_split/b???_* > /tmp/EXP00_namelist_rebuild
#    sdiff -w80 ./data/namelist /tmp/EXP00_namelist_rebuild
#
# Output contains lines which are before the first block, after the last
# block and lines between blocks.
#
# EVOLUTIONS
# ==========
#
# $Id$
#
# - fplod 2008-08-11T08:56:44Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#   * commentaires dans ce fichier en reStructuredText
#     cf. reStructuredText_ and Docutils_
#
# .. _reStructuredText: http://docutils.sourceforge.net/rst.html
# .. _Docutils: http://docutils.sourceforge.net/
#
# - fplod 2008-07-22T12:41:04Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#   * creation with test on
#     http://forge.ipsl.jussieu.fr/nemo/browser/trunk/CONFIG/ORCA2_LIM/EXP00/namelist revision 1151
#-
#
system=$(uname)
case "${system}" in
    AIX|IRIX64)
        echo " www : no specific posix checking"
    ;;
    *)
        #set -o posix
    ;;
esac
unset system
#
command=$(basename ${0})
log_date=$(date -u +"%Y%m%dT%H%M%SZ")
#
usage=" Usage : ${command} -i namelist -o dirout [-n]"
#
# default
blocknum=0
#
minargcount=4
if [ ${#} -lt ${minargcount} ]
then
    echo "eee : not enough arguments"
    echo "${usage}"
    exit 1
fi
#
while [ ${#} -gt 0 ]
do
    case ${1} in
        -i)
            filein=${2}
            shift
        ;;
        -o)
            dirout=${2}
            shift
        ;;
        -h)
            echo "${usage}"
            exit 0
        ;;
        -n)
            blocknum=1
        ;;
        *)
            # anything else
            echo "eee : unknown option ${1}"
            echo "eee : ${usage}"
            exit 1
        ;;
    esac
    # next flag
    shift
done
unset usage
#
set -u
#
# check for filein
if [ ! -r ${filein} ]
then
    echo "eee : ${filein} not accessible"
    exit 1
fi
#
# check for dirout
# the idea is to prevent unwanted overwrite
if [ -d ${dirout} ]
then
    echo "eee : ${dirout} already exist"
    exit 1
else
    mkdir -p ${dirout}
    if [ ${?} -ne 0 ]
    then
        echo "eee : cannot create ${dirout}"
        exit 1
    fi
fi
#
# loop on each line of ${filein} ie the namelist to be split
fileout=/dev/null
begin=0
if [ ${blocknum} -eq 1 ]
then
    iblock=0
fi
nbline=$(wc -l ${filein} | awk '{print $1}')
il=1
while [ ${il} -le ${nbline} ]
do
    line=$(sed -ne "${il},${il}p" ${filein})
    # if current line start with "&", we can establish the name of output file
    echo ${line} | grep -q "^&"
    if [ ${?} -eq 0 ]
    then
        begin=1
        # the form of this line is "&block   !"
        block=${line#&}
        block=${block%% *}
        echo "iii : detection of the beginning of ${block}"
        if [ ${blocknum} -eq 1 ]
        then
            iblock=$(( ${iblock} + 1))
            fileout=${dirout}/b$(echo ${iblock} | awk '{printf "%3.3d", $1}')_${block}
        else
            fileout="${dirout}/${block}"
        fi
        # check for ${fileout}
        if [ -f ${fileout} ]
        then
            echo "eee : ${fileout} already exists"
            exit 1
        fi
        # write the previous line and the current line
        echo "${memoline}" > ${fileout}
        echo "${line}" >> ${fileout}
    fi
    # if current line contains "/" there will be no more lines to put in fileout
    echo ${line} | grep -q "^/$"
    if [ ${?} -eq 0 ]
    then
        echo "iii : detection of the end of ${block}"
        echo "${line}" >> ${fileout}
        fileout="/dev/null"
    fi
    # if current line is between "&" and "/", it should be written in fileout
    # if current line is before the first &, it won't be written
    if [ "${fileout}" != "/dev/null" ]
    then
        if [ ${begin} -eq 0 ]
        then
            echo "${line}" >> ${fileout}
        fi
    fi
    # keep memory of this current line to use it if the next one is the
    # beginning of a block
    memoline=${line}
    begin=0
    # next line
    il=$(( ${il} + 1 ))
done
#
echo "iii : ${filein} is split in ${dirout}"
#
# end
exit 0
