#!/bin/sh

qsub -g tga-sca ./script/MASTER/reMASTER_optCDRG1.sh
qsub -g tga-sca ./script/MASTER/reMASTER_optCDRG2.sh
qsub -g tga-sca ./script/MASTER/reMASTER_optCDRG3.sh
qsub -g tga-sca ./script/MASTER/reMASTER_optCDRG4.sh
qsub -g tga-sca ./script/MASTER/reMASTER_optCDRG5.sh
