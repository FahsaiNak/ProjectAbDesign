#!/bin/sh

neo_resi=188
TCR="TCRb3"
[ -d "TCRgraft" ] || mkdir TCRgraft
[ -d "TCRgraft/${TCR}" ] || mkdir TCRgraft/${TCR}
for File in $(find ./TCRmatch/ -name "${TCR}_lowcutoff.match" -not -empty); do
	cdr=$(echo $File | cut -d "/" -f3)
	/gs/hs0/tga-sca/FAH/master-v1.6/bin/master --query ./PDS/${cdr}.pdb.pds --matchIn $File --structOut ./TCRmatch/${cdr}/${TCR}_lowcutoff.match.struct --outType match
	/gs/hs0/tga-sca/FAH/master-v1.6/bin/master --query ./PDS/${cdr}.pdb.pds --matchIn $File --structOut ./TCRmatch/${cdr}/${TCR}_lowcutoff.match.full.struct --outType full
	for matchpath in $(find ./TCRmatch/${cdr}/${TCR}_lowcutoff.match.struct/ -name "*.pdb"); do
		TCRpdb=$(head -1 $matchpath|rev| cut -d" "  -f2 | cut -d "/" -f1|rev| cut -d "_" -f1,2)
		pdbid=$(head -1 $matchpath|rev| cut -d" "  -f2 | cut -d "/" -f1|rev| cut -d "_" -f1)
		echo $TCRpdb $cdr $matchpath
		[ -d "./TCRgraft/${TCR}/${TCRpdb}" ] || mkdir ./TCRgraft/${TCR}/${TCRpdb}
		python ./script/CDRgraft.py "$matchpath" ../pdb/${cdr}.cdr.pdb $TCR ./TCRgraft/${TCR}/${TCRpdb}
		if [ ! -f "./TCRgraft/${TCR}/${TCRpdb}/${pdbid}_pMHC.pdb" ]; then
			cp /gs/hs0/tga-sca/FAH/human_tcr_STCRDab/imgt/${pdbid}.pdb ./TCRgraft/${TCR}/${TCRpdb}
			cp /gs/hs0/tga-sca/FAH/human_tcr_STCRDab/vb_fasta/${TCRpdb}_VB.fasta ./TCRgraft/${TCR}/${TCRpdb}
			python ./script/prepStruct.py /gs/hs0/tga-sca/FAH/rosettaMHC/p53R175H/rosettaMHC_5.pdb ./TCRgraft/${TCR}/${TCRpdb};
		fi;
		python ./script/TCRgraft_SASA.py $cdr $neo_resi ./TCRgraft/${TCR}/${TCRpdb};
	done;
done
printf "ID\tTCR\tShared\tNonshared\tNeositeContact\n" > ./TCRgraft/${TCR}/grafted_result.tsv # Nonshared = N lost interactions
for report in $(find ./TCRgraft -name interaction_report.txt); do
	cat $report >> ./TCRgraft/${TCR}/grafted_result.tsv;
done
	