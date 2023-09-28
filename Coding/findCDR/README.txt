# Step 0: Preparation #
## upload peptide-MHC structure in a single chain "A" to a working directory along with a "script" sub-folder ##
## create TCRmimic environment with MASTER program installed ##
## prepare PDB90 and AbAg databases along with all related files and change paths in the scripts ##
# Step 1: Epitope fragmentation and CDR-like motif identification #
> ./script/findCDRlike.sh
# Step 2: CDR-like candidate creation #
> ./script/CDRopt1.sh
> ./script/CDRopt2.sh
> ./script/CDRopt3.sh
> ./script/CDRopt4.sh
> find ./Ag* -name "*.cdr.pdb" > CDR-like_candidates_optimized.path
# Step 3: Solubility calculation
> mkdir solubility
> /gs/hs0/tga-sca/FAH/pdb2fasta/pdb2fasta2.sh $(find ./Ag* -name "*.cdr.pdb") > ./solubility/CDR-like_candidates.fasta
> cd solubility
## submit CDR-like_candidates.fasta to CamSol web server (https://www-cohsoftware.ch.cam.ac.uk/index.php) ##
## upload to solubility directory ##
> cat CamSol.txt| sed 1d| awk '$2 >= 1'| cut -f1| sort > CDR-like_candidates_highSol.name
> cat CamSol.txt| sed 1d| awk '$2 < 1'| cut -f1| sort > CDR-like_candidates_lowSol.name
> cat CamSol.txt| sed 1d| cut -f1| sort > CDR-like_candidates_withSol.name
> awk 'sub(/^>/, "")' CDR-like_candidates.fasta| sort > CDR-like_candidates.name
> comm -3 CDR-like_candidates.name CDR-like_candidates_withSol.name > CDR-like_candidates_noSol.name
> cp CDR-like_candidates_highSol.name CDR-like_candidates_Sol.name
> cat CDR-like_candidates_noSol.name >> CDR-like_candidates_Sol.name
> cat CDR-like_candidates_Sol.name | cut -d":" -f1 > CDR-like_candidates_Sol.cdr.name
> grep -f CDR-like_candidates_Sol.cdr.name ../CDR-like_candidates_optimized.path > CDR-like_candidates_Sol.path
> cd ..
# Step 4: neo-CDR-like candidate selection for CDR3 #
## select targeted antigen fragments; in this case, Ag_res5_185-189, Ag_res4_186-189 and Ag_res4_185-188 ##
> mkdir neoCDRcandidates
> grep "Ag_res5_185-189" ./solubility/CDR-like_candidates_Sol.path >> ./neoCDRcandidates/neoCDR-like_candidates_Sol.path
> grep "Ag_res4_186-189" ./solubility/CDR-like_candidates_Sol.path >> ./neoCDRcandidates/neoCDR-like_candidates_Sol.path
> grep "Ag_res4_185-188" ./solubility/CDR-like_candidates_Sol.path >> ./neoCDRcandidates/neoCDR-like_candidates_Sol.path
> mkdir pdb
> cat ./neoCDRcandidates/neoCDR-like_candidates_Sol.path | while read line;do cp $line ./neoCDRcandidates/pdb; done
> ./script/rankCDR.sh
## take a look at ranking.txt and select neo-CDR-like candidates ##
## selection criteria as follows ##
## 1) have an interaction with the neo-residue (True, NeositeContact) ##
## 2) high number of shared interactions (Shared) ##
## 3) low number of nonshared interactions (Nonshared) ##
