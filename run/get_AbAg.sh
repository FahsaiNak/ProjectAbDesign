#!/bin/bash
#set -e # stop on error
set -u # raise error if variable is unset
set -o pipefail # fail if any prior step failed

python ../src/process_AbAg.py --abag_filename ../Datasets/AbAg.pkl --output_filename ../Datasets/AbAg_processed.csv
