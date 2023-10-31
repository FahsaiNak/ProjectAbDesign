#!/bin/bash
#set -e # stop on error
set -u # raise error if variable is unset
set -o pipefail # fail if any prior step failed

get_info() {
    # Extract the PDB ID from the provided path.
    pdb=$(echo $1| cut -d "/" -f4| cut -d "." -f1)
    
    # Find match files containing the PDB ID and store the results in a temporary file.
    find ../Datasets/Ablike_top50 -name "*.match" | parallel grep -H ${pdb} > ${pdb}.match.tmp
    
    # Check if the temporary match file has content.
    if [ -s ${pdb}.match.tmp ]; then
        # Extract and sort unique CDR fragments from the match file.
        cat ${pdb}.match.tmp | cut -d":" -f1| sort -u > ${pdb}.cdr.tmp
        
        # Extract all matches and store in a temporary file.
        cat ${pdb}.match.tmp | cut -d":" -f2 > ${pdb}.all.match.tmp
        
        # Process CDR fragments one by one.
        cat ${pdb}.cdr.tmp| while read frag; do
            Name=$(echo $frag| cut -d "/" -f4| cut -d "." -f1)
            
            # Find common matches between the CDR fragment and all matches, then store in a temporary file.
            comm -12 <(sort ${pdb}.all.match.tmp) <(sort $frag) > ${frag}.${pdb}
            
            # Run a command (assuming it's a program) using the extracted CDR fragment and PDB ID.
            ../master-v1.6/bin/master --query ../Datasets/CDR_fragments_PDS/${Name}.pdb.pds --matchIn "${frag}.${pdb}" --structOut ../Datasets/AbAg/${Name}.match.${pdb}.struct --outType match --skipRMSD
            
            # Remove the temporary file.
            rm ${frag}.${pdb};
        done
        
        # Run a Python script to collect Ab-like information based on the PDB ID.
        python ../src/collect_Ablike.py --pdb $pdb --data_dir ../Datasets/AbAg/
        
        # Remove temporary structure files.
        rm -r ../Datasets/AbAg/*$pdb.struct;
    fi
    
    # Remove all temporary files related to the PDB ID.
    rm ${pdb}**.tmp
}
export -f get_info

# Find all PDB90 files with a specific extension and store in a temporary file.
find  ../Datasets/PDB90_PDS -type file -name "*.pdb.pds" > PDB90.pds.tmp

# Create the AbAg directory if it doesn't exist.
[ ! -d ../Datasets/AbAg ] && mkdir ../Datasets/AbAg

# Parallelize the process to get Ab-like information for multiple PDB files.
echo "getting Ab-like info"
parallel -j 2 -a PDB90.pds.tmp get_info

# Remove the temporary file containing PDB90 files.
rm PDB90.pds.tmp
