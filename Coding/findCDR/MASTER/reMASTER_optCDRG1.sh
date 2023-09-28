#!/bin/sh
#$-cwd
#$ -l q_core=1
#$ -l h_rt=24:00:00
#$ -N optCDR_G1

. /etc/profile.d/modules.sh
module load parallel
source ~/.bashrc
conda activate TCRmimic

doit() {
	Gpath="pathListaa"
	cdrepi=$(echo "$1" | awk -F/ '{print $NF}')
	qdir=$(echo "$1" | cut -d "/" -f1-4)
	echo $cdrepi
	head -50 ${qdir}/${Gpath}.match > ${qdir}/${Gpath}.top50.match
	/gs/hs0/tga-sca/FAH/master-v1.6/bin/master --query ${qdir}/*.pds --matchIn ${qdir}/${Gpath}.top50.match --structOut ${qdir}/${Gpath}.match.struct --outType match
	python ./script/optimizeCDR.py $qdir $Gpath ./6VR5_renumbered.pdb
	rm ${qdir}/${Gpath}.top50.match
	rm -r ${qdir}/${Gpath}.match.struct
}
export -f doit
parallel -j 10 -a CDR-epiContact.path doit