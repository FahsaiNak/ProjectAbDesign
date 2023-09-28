#!/bin/sh
#$-cwd
#$ -l q_core=1
#$ -l h_rt=24:00:00
#$ -N optCDR_G2

. /etc/profile.d/modules.sh
module load parallel
source ~/.bashrc
conda activate TCRmimic

doit() {
	Gpath="pathListab"
	TCR="TCRb3"
	contactpath=$(echo "$1")
	qdir=$(echo $contactpath|cut -d "/" -f1-6)
	pdb=$(echo $contactpath|cut -d "/" -f4)
	echo $contactpath $qdir $pdb
	/gs/hs0/tga-sca/FAH/master-v1.6/bin/master --query $contactpath.pds --matchIn ${qdir}/${Gpath}.match --structOut ${qdir}/${Gpath}.match.struct --outType match
	python ./script/optimizeTCR.py $qdir $Gpath ./TCRgraft/$TCR/$pdb/target_aligned.pdb
	rm -r ${qdir}/${Gpath}.match.struct
}
TCR="TCRb3"
export -f doit
parallel -j 10 -a ./TCRopt/${TCR}/TCR-epiContact.path doit