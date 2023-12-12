#!/bin/bash
#set -e # stop on error
set -u # raise error if variable is unset
set -o pipefail # fail if any prior step failed

parse_yaml() {
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s=\"%s\"\n", vn, $2, $3);
      }
   }'
}
export -f parse_yaml

findAglike() {
    # Print a message indicating the directory being searched.
    echo "searching in $1"
    
    # Run a Python script to find Ag-like regions using the provided pickle file.
    python ../src/find_Aglike.py --file $1
    
    # Remove the processed pickle file.
    rm $1
}
export -f findAglike

eval $(parse_yaml config.yaml)
# Find all pickle files in the AbAg directory and store them in a temporary file.
find $AbAg -type f -name "*.pkl" > PDB.pkl.tmp

# Print a message indicating the search for Ag-like regions.
echo "searching Ag-like regions"

# Parallelize the process to find Ag-like regions in multiple pickle files.
parallel -j $core findAglike :::: PDB.pkl.tmp

# Remove the temporary file containing the pickle file paths.
rm PDB.pkl.tmp

# Run a Python script to create a combined AbAg file from the processed data.
python ../src/create_AbAg.py --dir $AbAg

# Remove the entire AbAg directory.
rm -r $AbAg
