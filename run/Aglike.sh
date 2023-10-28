#!/bin/bash
#set -e # stop on error
set -u # raise error if variable is unset
set -o pipefail # fail if any prior step failed

findAglike() {
    echo "searching in $1"
    python ../src/find_Aglike.py --file $1
    rm $1
}
export -f findAglike

find  ../Datasets/AbAg -type file -name "*.pkl" > PDB.pkl.tmp
echo "searching Ag-like regions"
parallel -j 2 -a PDB.pkl.tmp findAglike
rm PDB.pkl.tmp
python ../src/create_AbAg.py --dir ../Datasets/AbAg
rm -r ../Datasets/AbAg
