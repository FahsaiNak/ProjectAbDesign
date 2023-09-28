#!/bin/sh

qsub -g tga-sca ./script/MASTER/MASTER_optTCRG1.sh
qsub -g tga-sca ./script/MASTER/MASTER_optTCRG2.sh
qsub -g tga-sca ./script/MASTER/MASTER_optTCRG3.sh
qsub -g tga-sca ./script/MASTER/MASTER_optTCRG4.sh
qsub -g tga-sca ./script/MASTER/MASTER_optTCRG5.sh
