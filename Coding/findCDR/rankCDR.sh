#!/bin/sh
#$-cwd
#$ -l q_core=1
#$ -l h_rt=01:00:00
#$ -N rankOptCDR

. /etc/profile.d/modules.sh
module load parallel
source ~/.bashrc
conda activate TCRmimic

File="./neoCDRcandidates/neoCDR-like_candidates_Sol.path"
doit() {
	neo_resi=188
	epi=$(echo "$1" |rev | cut -d "/" -f3 | rev)
	code_id=$(echo "$1" | awk -F/ '{print $NF}' | cut -d "." -f1)
	abag_id=$(echo "$code_id" | cut -d "_" -f1)
	ag_pos1=$(echo "$code_id" | cut -d"_" -f2 | awk -F 'Ag' '{ print $2}' | cut -d "-" -f1)
	ag_pos2=$(echo "$code_id" | cut -d"_" -f2 | awk -F 'Ag' '{ print $2}' | cut -d "-" -f2)
	pdb_id=$(grep "$abag_id" CDR-like_candidates.lst| grep "$epi"| grep "${ag_pos1},"| grep " ${ag_pos2}"| cut -f2)
	cdr_chain=$(grep "$abag_id" CDR-like_candidates.lst| grep "$epi"| grep "${ag_pos1},"| grep " ${ag_pos2}"| cut -f3)
	cdr_list=$(grep "$abag_id" CDR-like_candidates.lst| grep "$epi"| grep "${ag_pos1},"| grep " ${ag_pos2}"| cut -f4)
	ag_chain=$(grep "$abag_id" CDR-like_candidates.lst| grep "$epi"| grep "${ag_pos1},"| grep " ${ag_pos2}"| cut -f5)
	ag_list=$(grep "$abag_id" CDR-like_candidates.lst| grep "$epi"| grep "${ag_pos1},"| grep " ${ag_pos2}"| cut -f6)
	epifrag_list=$(grep "$abag_id" CDR-like_candidates.lst| grep "$epi"| grep "${ag_pos1},"| grep " ${ag_pos2}"| cut -f8)
	echo $code_id
	qdir="./neoCDRcandidates"
	python ./script/rankCDR.py $code_id /gs/hs0/tga-sca/FAH/all_PDB/ncAA/${pdb_id}.pdb $cdr_chain "$cdr_list" $ag_chain "$ag_list" $epi "$epifrag_list" ./6VR5_renumbered.pdb $neo_resi $qdir
}
printf "ID\tEpitope\tShared\tNonshared\tNeositeContact\n" > ./neoCDRcandidates/ranking.txt
export -f doit
parallel -j 10 -a $File doit