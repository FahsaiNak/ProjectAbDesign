#!/bin/sh
#$-cwd
#$ -l q_core=1
#$ -l h_rt=00:10:00
#$ -N optCDR

. /etc/profile.d/modules.sh
module load parallel
source ~/.bashrc
conda activate TCRmimic

doit() {
	python ./script/combindOptCDR.py $1
}
export -f doit
find ./Ag* -type d -name "*_EpiAb" -not -empty > optimized_CDR.dir
parallel -j 10 -a optimized_CDR.dir doit