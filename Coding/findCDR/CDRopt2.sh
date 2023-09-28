#!/bin/sh

qsub -g tga-sca ./script/MASTER/MASTER_optCDRG1.sh
qsub -g tga-sca ./script/MASTER/MASTER_optCDRG2.sh
qsub -g tga-sca ./script/MASTER/MASTER_optCDRG3.sh
qsub -g tga-sca ./script/MASTER/MASTER_optCDRG4.sh
qsub -g tga-sca ./script/MASTER/MASTER_optCDRG5.sh
