#!/bin/bash
set -e # stop on error
set -u # raise error if variable is unset
set -o pipefail # fail if any prior step failed

extract_CDR_from_PDB_file=extract_CDRs_from_PDB.py
fragment_CDR_file=fragment_CDRs.py
python ../src/$extract_CDR_from_PDB_file
python ../src/$fragment_CDR_file
