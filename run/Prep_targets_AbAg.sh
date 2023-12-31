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

prepTargets() {
   pdb=$(echo $1| cut -d "/" -f4)
	$2/createPDS --type target --pdb $1 --pds $3/$pdb.pds
}
export -f prepTargets

eval $(parse_yaml config.yaml)
find $PDB90 -type f -name "*.pdb" > PDB.lst.tmp
echo "The targets are processing..."
[ ! -d $PDS90 ] && mkdir $PDS90
parallel -j $core prepTargets :::: PDB.lst.tmp ::: $MASTER ::: $PDS90
rm PDB.lst.tmp
