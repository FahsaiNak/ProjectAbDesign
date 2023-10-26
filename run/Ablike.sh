#!/bin/bash
#set -e # stop on error
#set -u # raise error if variable is unset
#set -o pipefail # fail if any prior step failed

findAblike() {
    echo "$1 start!"
    Name=$(echo $1| cut -d "/" -f4| cut -d "." -f1)
    Numres=$(grep -a " CA " $1 | cat -n | tail -1| awk -F ' ' '{ print $1}')
    V=$(echo "scale=2;(${Numres}*0.05)+0.4" | bc)
    if [[ $V > 1.00 ]]; then V=1.00; fi
    ../master-v1.6/bin/master --query $1 --targetList PDB90.pds.lst.tmp --rmsdCut $V --matchOut ../Datasets/Ablike/${Name}.match
    head -50 ../Datasets/Ablike/${Name}.match > ../Datasets/Ablike_top50/${Name}.match
    rm ../Datasets/Ablike/${Name}.match
}
export -f findAblike

find  ../Datasets/PDB90_PDS -type file -name "*.pdb.pds" > PDB90.pds.lst.tmp
find  ../Datasets/CDR_fragments_PDS -type file -name "*.pdb.pds"| sort > CDR.frag.pds.lst.tmp
echo "Searching for CDR-like regions ..."
[ ! -d ../Datasets/Ablike ] && mkdir ../Datasets/Ablike
[ ! -d ../Datasets/Ablike_top50 ] && mkdir ../Datasets/Ablike_top50
parallel -j 3 -a CDR.frag.pds.lst.tmp findAblike
rm *lst.tmp
rm -r ../Datasets/Ablike
find ../Datasets/Ablike_top50/ -name "*.match" -empty -delete
