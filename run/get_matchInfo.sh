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

get_info() {
    # Extract the PDB ID from the provided path.
    pdb=$(echo $1| cut -d "/" -f4| cut -d "." -f1)
    
    # Find match files containing the PDB ID and store the results in a temporary file.
    find $5 -name "*.match" | parallel grep -H ${pdb} > ${pdb}.match.tmp
    
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
            $2/master --query $4/${Name}.pdb.pds --matchIn "${frag}.${pdb}" --structOut $3/${Name}.match.${pdb}.struct --outType match --skipRMSD
            
            # Remove the temporary file.
            rm ${frag}.${pdb};
        done
        
        # Run a Python script to collect Ab-like information based on the PDB ID.
        python ../src/collect_Ablike.py --pdb $pdb --data_dir $3
        
        # Remove temporary structure files.
        rm -r $3/**$pdb.struct;
    fi
    
    # Remove all temporary files related to the PDB ID.
    rm ${pdb}**.tmp
}
export -f get_info

eval $(parse_yaml config.yaml)
# Find all PDB90 files with a specific extension and store in a temporary file.
find $PDS90 -type f -name "*.pdb.pds" > PDB90.pds.tmp

# Create the AbAg directory if it doesn't exist.
[ ! -d $AbAg ] && mkdir $AbAg

# Parallelize the process to get Ab-like information for multiple PDB files.
echo "getting Ab-like info"
parallel -j $core get_info :::: PDB90.pds.tmp ::: $MASTER ::: $AbAg ::: $PDSCDRfrag ::: $Ablike

# Remove the temporary file containing PDB90 files.
rm PDB90.pds.tmp
