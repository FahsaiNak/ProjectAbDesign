#!/bin/sh
#$-cwd
#$ -l q_core=1
#$ -l h_rt=24:00:00
#$ -N findCDRlike

. /etc/profile.d/modules.sh
module load parallel
source ~/.bashrc
conda activate TCRmimic

python ./script/linearFrag.py ./6VR5_renumbered.pdb ./
find -type d -name "Ag*" | cut -d "/" -f2 > epitope_frag.name
echo "Fragmentation Done"

doit() {
	Name=$(echo $1)
	echo $Name
	/gs/hs0/tga-sca/FAH/pdb2fasta/pdb2fasta2.sh ./${Name}/${Name}.pdb > ./${Name}/${Name}.fasta
	blastp -query ./${Name}/${Name}.fasta -db /gs/hs0/tga-sca/FAH/AbAg/p53R175H/Ag-like/fasta/Ag-like.fa -qcov_hsp_perc 100.0 -matrix BLOSUM62 -task 'blastp-short' -word_size 2 -seg 'no' -evalue 20000 -ungapped -comp_based_stats F -max_target_seqs 60000 -outfmt 6 -out ./${Name}/blast_hits.txt
    /gs/hs0/tga-sca/FAH/master-v1.6/bin/createPDS --type query --pdb ./${Name}/${Name}.pdb --pds ./${Name}/${Name}.pdb.pds
    cat ./${Name}/blast_hits.txt | cut -f2 | cut -d ":" -f1 > ./${Name}/Ag-like.name
    echo "Sequence Matching Done"
    Numres=$(grep -a " CA " ./${Name}/${Name}.pdb | cat -n | tail -1| awk -F ' ' '{ print $1}')
    V=$(echo "scale=2;(${Numres}*0.033)+0.4" | bc)
	grep -f ./${Name}/Ag-like.name /gs/hs0/tga-sca/FAH/AbAg/p53R175H/Ag-like/p53R175H_linear.match.targetList > ./${Name}/TargetList.path
	/gs/hs0/tga-sca/FAH/master-v1.6/bin/master --query ./${Name}/${Name}.pdb.pds --targetList ./${Name}/TargetList.path --rmsdCut $V --matchOut ./${Name}/${Name}.match --structOut ./${Name}/${Name}.match.struct --outType match;
	echo "Structure Matching Done"
	if [ -s ./${Name}/${Name}.match ]; then
		touch ./${Name}/CDR-like_candidates.lst
		for match in $(find ./${Name}/${Name}.match.struct -type f -name "*.pdb");do
			Ag=$(head -1 $match| awk -F ' ' '{ print $3}'| awk -F/ '{print $NF}' | cut -d "." -f1)
			pdb=$(echo $Ag| cut -d "-" -f1)
			chain=$(echo $Ag| cut -d "-" -f2| cut -d "_" -f1)
			pos1=$(echo $Ag| cut -d "-" -f2| cut -d "_" -f2)
			pos2=$(echo $Ag| cut -d "-" -f3)
			echo $Ag $pdb $chain $pos1 $pos2
			grep $pdb /gs/hs0/tga-sca/FAH/AbAg/p53R175H/CDR-Ag.lst | grep "$chain" | grep "${pos1}, ${pos2}]]"| sort -u >> ./${Name}/CDR-Ag.match
			grep $pdb ./${Name}/CDR-Ag.match | grep "$chain"| grep "${pos1}, ${pos2}]]" > ./${Name}/CDR-Ag.one.match
			awk 'FNR==NR { a[$1]=$1; next } $1 in a { print $0}' ./${Name}/CDR-Ag.one.match /gs/hs0/tga-sca/FAH/AbAg/p53R175H/AbAg.lst > ./${Name}/CDR-Ag.one.lst
			python ./script/findClashes.py $Name $match ./${Name}/CDR-Ag.one.lst ./6VR5_renumbered.pdb ./${Name}
			rm ./${Name}/CDR-Ag.one.*
			rm $match;
		done;
	fi
	rm -r ./${Name}/${Name}.match.struct
}
export -f doit
parallel -j 10 -a epitope_frag.name doit

for path in $(find ./Ag_* -name "CDR-Ag.match" -not -empty);do cat $path >> CDR-Ag.match;done
sort -u CDR-Ag.match > CDR-Ag.match.tmp
mv CDR-Ag.match.tmp CDR-Ag.match
for path in $(find ./Ag_* -name "CDR-like_candidates.lst" -not -empty);do cat $path >> CDR-like_candidates.lst;done
sort -u CDR-like_candidates.lst > CDR-like_candidates.lst.tmp
mv CDR-like_candidates.lst.tmp CDR-like_candidates.lst
echo "CDR-like Identification Done"


