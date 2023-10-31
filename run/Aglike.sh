#!/bin/bash
#set -e # stop on error
set -u # raise error if variable is unset
set -o pipefail # fail if any prior step failed

findAglike() {
    # Print a message indicating the directory being searched.
    echo "searching in $1"
    
    # Run a Python script to find Ag-like regions using the provided pickle file.
    python ../src/find_Aglike.py --file $1
    
    # Remove the processed pickle file.
    rm $1
}
export -f findAglike

# Find all pickle files in the AbAg directory and store them in a temporary file.
find  ../Datasets/AbAg -type file -name "*.pkl" > PDB.pkl.tmp

# Print a message indicating the search for Ag-like regions.
echo "searching Ag-like regions"

# Parallelize the process to find Ag-like regions in multiple pickle files.
parallel -j 2 -a PDB.pkl.tmp findAglike

# Remove the temporary file containing the pickle file paths.
rm PDB.pkl.tmp

# Run a Python script to create a combined AbAg file from the processed data.
python ../src/create_AbAg.py --dir ../Datasets/AbAg

# Remove the entire AbAg directory.
rm -r ../Datasets/AbAg
