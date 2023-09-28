#!/bin/sh
#$-cwd
#$ -l q_core=1
#$ -l h_rt=03:00:00
#$ -N optCDR

. /etc/profile.d/modules.sh
module load parallel
source ~/.bashrc
conda activate TCRmimic

doit() {
	neo_resi=188
	tcrpath=$1
	qdir=$(echo $tcrpath|cut -d "/" -f1-5)
	opt_tcr_code=$(echo $tcrpath|rev|cut -d "/" -f1|rev|cut -d "." -f1)
	python ./script/rankOptTCR.py $qdir $opt_tcr_code $neo_resi
}
export -f doit

TCR="TCRb3"
cat ./TCRopt/${TCR}/TCR-epiContact.path|cut -d "/" -f1-5|sort -u > ./TCRopt/${TCR}/TCR-opt.path
cat ./TCRopt/${TCR}/TCR-opt.path| while read line; do
	python ./script/combindOptTCR.py $line
	cp TCRgraft/${TCR}/$(echo $line|cut -d "/" -f4)/*_pMHC.pdb $(echo $line|cut -d "/" -f1-4)
	ls $line/*.pdb > ./TCRopt/${TCR}/optimizedTCRm_candidates.name
	parallel -j 10 -a ./TCRopt/${TCR}/optimizedTCRm_candidates.name doit
	rm ./TCRopt/${TCR}/optimizedTCRm_candidates.name
done
