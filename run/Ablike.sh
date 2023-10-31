#!/bin/bash
#set -e # stop on error
set -u # raise error if variable is unset
set -o pipefail # fail if any prior step failed

findAblike() {
    # Print a message indicating the start of processing for the provided file.
    echo "$1 start!"
    
    # Extract the name and the number of residues from the provided file.
    Name=$(echo $1| cut -d "/" -f4| cut -d "." -f1)
    Numres=$(grep -a " CA " $1 | cat -n | tail -1| awk -F ' ' '{ print $1}')
    
    # Calculate a value 'V' based on the number of residues.
    V=$(echo "scale=2;(${Numres}*0.05)+0.4" | bc)
    
    # Ensure that 'V' is not greater than 1.00.
    if [[ $V > 1.00 ]]; then
        V=1.00
    fi
    
    # Run a command (assuming it's a program) using the provided file and calculated 'V'.
    ../master-v1.6/bin/master --query $1 --targetList PDB90.pds.lst.tmp --rmsdCut $V --matchOut ../Datasets/Ablike/${Name}.match
    
    # Extract the top 50 lines from the generated match file and save it to a separate file.
    head -50 ../Datasets/Ablike/${Name}.match > ../Datasets/Ablike_top50/${Name}.match
    
    # Remove the original match file.
    rm ../Datasets/Ablike/${Name}.match
}
export -f findAblike

# Find all PDB90 files with a specific extension and store them in a temporary file.
find  ../Datasets/PDB90_PDS -type file -name "*.pdb.pds" > PDB90.pds.lst.tmp

# Find all CDR fragment files and sort them.
find  ../Datasets/CDR_fragments_PDS -type file -name "*.pdb.pds"| sort > CDR.frag.pds.lst.tmp

# Print a message indicating the search for CDR-like regions.
echo "Searching for CDR-like regions ..."

# Create Ablike and Ablike_top50 directories if they don't exist.
[ ! -d ../Datasets/Ablike ] && mkdir ../Datasets/Ablike
[ ! -d ../Datasets/Ablike_top50 ] && mkdir ../Datasets/Ablike_top50

# Parallelize the process to find CDR-like regions in multiple files.
parallel -j 3 -a CDR.frag.pds.lst.tmp findAblike

# Remove temporary files containing file lists.
rm *lst.tmp

# Remove the entire Ablike directory.
rm -r ../Datasets/Ablike

# Delete empty match files in the Ablike_top50 directory.
find ../Datasets/Ablike_top50/ -name "*.match" -empty -delete
