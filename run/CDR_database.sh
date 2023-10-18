#!/bin/bash
set -e # stop on error
set -u # raise error if variable is unset
set -o pipefail # fail if any prior step failed

python ../src/CDR_frag.py --path ../Datasets/CDR --savePath ../Datasets/CDR_fragments
