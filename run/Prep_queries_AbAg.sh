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

prepQueries() {
	pdb=$(echo $1| cut -d "/" -f4)
	$2/createPDS --type query --pdb $1 --pds $3/$pdb.pds
}
export -f prepQueries

eval $(parse_yaml config.yaml)
find  $PDBCDRfrag -type file -name "*.pdb" > CDR.lst.tmp
echo "The queries are processing..."
[ ! -d $PDSCDRfrag ] && mkdir $PDSCDRfrag
parallel -j $core prepQueries :::: CDR.lst.tmp ::: $MASTER ::: $PDSCDRfrag
rm CDR.lst.tmp
