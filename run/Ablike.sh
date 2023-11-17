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

findAblike() {
    # Print a message indicating the start of processing for the provided file.
    echo "$1 start!"
    
    # Extract the name and the number of residues from the provided file.
    Name=$(echo $1| cut -d "/" -f4| cut -d "." -f1)
    Numres=$(grep -a " CA " $1 | cat -n | tail -1| awk -F ' ' '{ print $1}')
    
    # Calculate a value 'V' based on the number of residues.
    V=$(echo "scale=2;(${Numres}*0.05)+0.4" | bc)
    
    # Ensure that 'V' is not greater than 1.00.
    if [[ $V > $5 ]]; then
        V=$5
    fi
    
    # Run a command (assuming it's a program) using the provided file and calculated 'V'.
    $2/master --query $1 --targetList PDB90.pds.lst.tmp --rmsdCut $V --matchOut $3/${Name}.match
    
    # Extract the top 50 lines from the generated match file and save it to a separate file.
    head -${4} $3/${Name}.match > "${3}_top${4}/${Name}.match"
    
    # Remove the original match file.
    rm $3/${Name}.match
}
export -f findAblike
eval $(parse_yaml config.yaml)

# Find all PDB90 files with a specific extension and store them in a temporary file.
find $PDS90 -type file -name "*.pdb.pds" > PDB90.pds.lst.tmp

# Find all CDR fragment files and sort them.
find $PDSCDRfrag -type file -name "*.pdb.pds"| sort > CDR.frag.pds.lst.tmp

# Print a message indicating the search for CDR-like regions.
echo "Searching for CDR-like regions ..."

# Create Ablike and Ablike_top50 directories if they don't exist.
[ ! -d $Ablike ] && mkdir $Ablike
[ ! -d ${Ablike}_top${topselect} ] && mkdir ${Ablike}_top${topselect}

# Parallelize the process to find CDR-like regions in multiple files.
parallel -j $core findAblike :::: CDR.frag.pds.lst.tmp ::: $MASTER ::: $Ablike ::: $topselect ::: $maxcutoff

# Remove temporary files containing file lists.
rm *lst.tmp

# Remove the entire Ablike directory.
rm -r $Ablike

# Delete empty match files in the Ablike_top50 directory.
find "${Ablike}_top${topselect}" -name "*.match" -empty -delete
