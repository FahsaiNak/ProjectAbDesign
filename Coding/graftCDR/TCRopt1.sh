#!/bin/sh
#$-cwd
#$ -l q_core=1
#$ -l h_rt=24:00:00
#$ -N TCRopt1

. /etc/profile.d/modules.sh
module load parallel
source ~/.bashrc
conda activate TCRmimic

TCR="TCRb3"
[ -d "./TCRopt" ] || mkdir TCRopt
[ -d "./TCRopt/${TCR}" ] || mkdir TCRopt/${TCR}
grep "True" $(find ./TCRgraft/${TCR} -name "grafted_result.tsv") > ./TCRopt/${TCR}/TCRm_candidate.lst
doit() {
	echo $1
	TCR="TCRb3"
	pdb=$(echo $1 | cut -d " " -f2)
	tcr_code=$(echo $1 | cut -d " " -f1)
	[ -d "./TCRopt/${TCR}/${pdb}" ] || mkdir TCRopt/${TCR}/${pdb}
	[ -d "./TCRopt/${TCR}/${pdb}/${tcr_code}" ] || mkdir TCRopt/${TCR}/${pdb}/${tcr_code}
	python ./script/Epicontact_SASA.py "$1" $TCR
}
export -f doit
parallel -a TCRopt/${TCR}/TCRm_candidate.lst doit
find ./TCRopt/${TCR} -name "*contact.pdb" > ./TCRopt/${TCR}/TCR-epiContact.path


