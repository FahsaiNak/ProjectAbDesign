#!/bin/bash
#set -e # stop on error
#set -u # raise error if variable is unset
#set -o pipefail # fail if any prior step failed

get_info() {
    pdb=$(echo $1| cut -d "/" -f4| cut -d "." -f1)
    find ../Datasets/Ablike_top50_old -name "*.match" | parallel grep -H ${pdb} > ${pdb}.match.tmp
    if [ -s ${pdb}.match.tmp ]; then
        cat ${pdb}.match.tmp | cut -d":" -f1| sort -u > ${pdb}.cdr.tmp
        cat ${pdb}.match.tmp | cut -d":" -f2 > ${pdb}.all.match.tmp
        cat ${pdb}.cdr.tmp| while read frag; do
            Name=$(echo $frag| cut -d "/" -f4| cut -d "." -f1)
            comm -12 <(sort ${pdb}.all.match.tmp) <(sort $frag) > ${frag}.${pdb}
            ../master-v1.6/bin/master --query ../Datasets/CDR_fragments_PDS/${Name}.pdb.pds --matchIn "${frag}.${pdb}" --structOut ../Datasets/AbAg/${Name}.match.${pdb}.struct --outType match --skipRMSD
            rm ${frag}.${pdb};
        done
        python ../src/collect_Ablike.py --pdb $pdb --data_dir ../Datasets/AbAg/;
    fi
    rm ${pdb}**.tmp
}
export -f get_info

find  ../Datasets/PDB90_PDS -type file -name "*.pdb.pds" > PDB90.pds.tmp
[ ! -d ../Datasets/AbAg ] && mkdir ../Datasets/AbAg
echo "getting Ab-like info"
parallel -j 2 -a PDB90.pds.tmp get_info
rm PDB90.pds.tmp
