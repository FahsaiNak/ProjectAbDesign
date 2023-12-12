#!/bin/bash
#set -e # stop on error
set -u # raise error if variable is unset
set -o pipefail # fail if any prior step failed

python ../src/PDB90.py --output_folder "../Datasets/all_PDB" --csv_file "../Datasets/somePDB.csv"