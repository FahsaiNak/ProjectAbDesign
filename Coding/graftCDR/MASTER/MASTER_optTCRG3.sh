#!/bin/sh
#$-cwd
#$ -l q_core=1
#$ -l h_rt=24:00:00
#$ -N pl_MASTER_G3

. /etc/profile.d/modules.sh
module load parallel

doit() {
	Gpath="pathListac"
	contactpath=$(echo "$1")
	qdir=$(echo "$contactpath"|cut -d "/" -f2-6)
	[ -s $contactpath.pds ] || /gs/hs0/tga-sca/FAH/master-v1.6/bin/createPDS --type query --pdb $contactpath --pds $contactpath.pds
	Numres=$(grep -a " CA " $contactpath.pds | cat -n | tail -1| awk -F ' ' '{ print $1}')
	V=$(echo "scale=2;(${Numres}*0.025)+0.6" | bc)
	echo $contactpath $V
	/gs/hs0/tga-sca/FAH/master-v1.6/bin/master --query $contactpath.pds --targetList /gs/hs0/tga-sca/FAH/all_PDB/split/10K/${Gpath} --rmsdCut $V --matchOut ${qdir}/${Gpath}.match
}
TCR="TCRb3"
export -f doit
parallel -j 10 -a ./TCRopt/${TCR}/TCR-epiContact.path doit
