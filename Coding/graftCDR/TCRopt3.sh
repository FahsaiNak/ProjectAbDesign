#!/bin/sh

qsub -g tga-sca ./script/MASTER/reMASTER_optTCRG1.sh
qsub -g tga-sca ./script/MASTER/reMASTER_optTCRG2.sh
qsub -g tga-sca ./script/MASTER/reMASTER_optTCRG3.sh
qsub -g tga-sca ./script/MASTER/reMASTER_optTCRG4.sh
qsub -g tga-sca ./script/MASTER/reMASTER_optTCRG5.sh
