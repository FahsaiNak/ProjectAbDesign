#!/bin/sh
#$-cwd
#$ -l q_core=1
#$ -l h_rt=24:00:00
#$ -N CDRopt1

. /etc/profile.d/modules.sh
module load parallel
source ~/.bashrc
conda activate TCRmimic

doit() {
	line=$(echo $1)
	echo $line
	python ./script/Epicontact_SASA.py "$line" 
}
export -f doit
parallel -j 10 -a CDR-like_candidates.lst doit
find ./Ag_* -name "*contact.pdb" > ./CDR-epiContact.path


