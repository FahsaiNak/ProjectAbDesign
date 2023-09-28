#!/bin/sh
#$-cwd
#$ -l q_core=1
#$ -l h_rt=24:00:00
#$ -N pl_MASTER_G2

. /etc/profile.d/modules.sh
module load parallel

doit() {
	Gpath="pathListab"
	epi=$(echo "$1" | cut -d "/" -f2)
	cdr=$(echo "$1" | cut -d "/" -f3)
	cdrepi=$(echo "$1" | cut -d "/" -f4)
	qdir="${epi}/${cdr}/${cdrepi}"
	[ -s ${qdir}/${cdrepi}.pdb.pds ] || /gs/hs0/tga-sca/FAH/master-v1.6/bin/createPDS --type query --pdb ${qdir}/${cdrepi}.pdb --pds ${qdir}/${cdrepi}.pdb.pds
	Numres=$(grep -a " CA " ${qdir}/${cdrepi}.pdb.pds | cat -n | tail -1| awk -F ' ' '{ print $1}')
	V=$(echo "scale=2;(${Numres}*0.025)+0.6" | bc)
	echo $epi $cdrepi $V
	/gs/hs0/tga-sca/FAH/master-v1.6/bin/master --query ${qdir}/${cdrepi}.pdb.pds --targetList /gs/hs0/tga-sca/FAH/all_PDB/split/10K/${Gpath} --rmsdCut $V --matchOut ${qdir}/${Gpath}.match
}

export -f doit
parallel -j 10 -a CDR-epiContact.lst doit
