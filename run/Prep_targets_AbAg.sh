#!/bin/bash
#set -e # stop on error
set -u # raise error if variable is unset
set -o pipefail # fail if any prior step failed

prepTargets() {
	pdb=$(echo $1| cut -d "/" -f4)
	../master-v1.6/bin/createPDS --type target --pdb $1 --pds ../Datasets/PDB90_PDS/$pdb.pds
}
export -f prepTargets

find  ../Datasets/all_PDB -type file -name "*.pdb" > PDB.lst.tmp
echo "The targets are processing..."
[ ! -d ../Datasets/all_PDB_PDS ] && mkdir ../Datasets/PDB90_PDS
parallel -j 5 -a PDB.lst.tmp prepTargets
rm PDB.lst.tmp
