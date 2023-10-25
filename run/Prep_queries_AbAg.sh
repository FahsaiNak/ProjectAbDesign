#!/bin/bash
#set -e # stop on error
set -u # raise error if variable is unset
set -o pipefail # fail if any prior step failed

prepQueries() {
	pdb=$(echo $1| cut -d "/" -f4)
	../master-v1.6/bin/createPDS --type query --pdb $1 --pds ../Datasets/CDR_fragments_PDS/$pdb.pds
}
export -f prepQueries

find  ../Datasets/CDR_fragments -type file -name "*.pdb" > CDR.lst.tmp
echo "The queries are processing..."
[ ! -d ../Datasets/CDR_fragments_PDS ] && mkdir ../Datasets/CDR_fragments_PDS
parallel -j 5 -a CDR.lst.tmp prepQueries
rm CDR.lst.tmp
