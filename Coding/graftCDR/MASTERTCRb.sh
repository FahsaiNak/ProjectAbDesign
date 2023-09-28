[ -d PDS ] || mkdir PDS
[ -d TCRmatch ] || mkdir TCRmatch
cat shortlist.cdr.name | while read Name; do
	[ -d ./TCRmatch/${Name} ] || mkdir ./TCRmatch/${Name}
	/gs/hs0/tga-sca/FAH/master-v1.6/bin/createPDS --type query --pdb ../pdb/${Name}.cdr.pdb --pds ./PDS/${Name}.pdb.pds
	Numres=$(grep -a " CA " ./PDS/${Name}.pdb.pds| cat -n | tail -1| awk -F ' ' '{ print $1}')
	V=$(echo "scale=2;(${Numres}*0.025)+0.6" | bc) #scale=2;(${Numres}*0.05)+0.4
	echo $Name $Numres $V
	chain="b3"
	/gs/hs0/tga-sca/FAH/master-v1.6/bin/master --query ./PDS/${Name}.pdb.pds --targetList /gs/hs0/tga-sca/FAH/human_tcr_STCRDab/PDS/CDR_${chain}.path --rmsdCut $V --bbRMSD --matchOut ./TCRmatch/${Name}/TCR${chain}_lowcutoff.match --bbRMSD
	chain="b1"
	/gs/hs0/tga-sca/FAH/master-v1.6/bin/master --query ./PDS/${Name}.pdb.pds --targetList /gs/hs0/tga-sca/FAH/human_tcr_STCRDab/PDS/CDR_${chain}.path --rmsdCut $V --bbRMSD --matchOut ./TCRmatch/${Name}/TCR${chain}_lowcutoff.match --bbRMSD;
	find ./TCRmatch/${Name}/*_lowcutoff.match -empty -delete;
done
