# Step 0: Preparation #
## create subdirectory "shortlist" in "neoCDRcandidates" ##
## upload a text file of selected neo-CDR-like candidate names ##
## TCR structure (.pdb) dataset including fasta file of TCR variable beta chain and beta-CDR3 fragments (.pdb and .pds) of each ##
## find here: http://opig.stats.ox.ac.uk/webapps/stcrdab/ ##
## activate environment ##
# Step 1: CDR grafting #
> ./script/MASTERTCRb3.sh
> ./script/TCRgraft.sh
# Step 2: CDR-grafted-TCRm optimization #
> ./script/TCRopt1.sh #find contacting-antigen residues for each CDR-grafted-TCRvb3 residue
> ./script/TCRopt2.sh #search for MASTER-matching structures of contacting-antigen-TCRm from PDB90
> ./script/TCRopt3.sh #replace the MASTER-matching residues to the original CDR-grafted TCR
> qsub -g tga-sca ./script/TCRopt4.sh #combine all potential residues

