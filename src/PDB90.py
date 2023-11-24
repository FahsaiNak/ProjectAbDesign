"""
This script is used for downloading, renaming,
and cleaning protein data bank (PDB) files.

Modules and Packages:
- os: Provides functions for interacting with the operating system.
- time: Provides various time-related functions.
- glob: Finds all the pathnames matching a specified pattern.
- subprocess: Allows you to spawn new processes, connect to their
input/output/error pipes, and obtain their return codes.
- argparse: Makes it easy to write user-friendly command-line interfaces.
- Bio.PDB: A Biopython module that focuses on working with crystal
structures of biological macromolecules.

Classes:
- NotHetero: A class that inherits from the Select class in the Bio.PDB module.
It is used to select only the non-hetero parts of a protein by removing
small molecules or ions that are not a main part of the protein.

Functions:
- download_and_uncompress_files(folder, csv_file): Downloads and uncompresses
PDB files from a given CSV file into a specified folder.
- rename_files(directory): Renames all files in a specified directory with
a .ent termination to .pdb.
- clean_all_pdb_files(folder): Cleans all PDB files in a specified
folder by keeping only the non-hetero parts of a protein.

Main Function:
The main() function is the entry point of the script. It parses command-line
arguments for the output folder and CSV file, and then calls the
download_and_uncompress_files(), rename_files(),
and clean_all_pdb_files() functions.

Execution:
The script is executed from the command line with the following arguments:
- --output_folder: The directory to save the cleaned files.
- --csv_file: The CSV file to use.

Example usage:
python script.py --output_folder /path/to/output/folder
--csv_file /path/to/csv/file.csv

Please replace /path/to/output/folder and
/path/to/csv/file.csv with your actual paths.

This script is designed to be used with Python 3.
Please ensure that you have the necessary
Python packages installed before running the script.
You can install the necessary packages
using pip or by creating an environment as stated in the README.md file.

"""

import os
import glob
import subprocess
import argparse
from Bio.PDB import PDBParser, PDBIO, Select


# Canonical Aminoacid residues to keep in PDB sequences
CANONICAL_RESIDUES = [
    'ALA', 'ARG', 'ASN', 'ASP',
    'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS',
    'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL'
]


class NotHetero(Select):
    # This function keeps only the non-hetero parts of a protein
    # by removing small molecules or ions
    # that are not a main part of the protein
    def accept_residue(self, residue):
        # If the sequence contains water molecules (HOH or WAT)
        # return 0, removing them
        if residue.get_resname() in ['HOH', 'WAT']:
            return 0
        # Keep only 1-meres, 2-meres, and 3-meres in the sequence
        elif len(residue.get_resname()) > 3:
            return 0
        # If the residue 3-mere is a canonical AA residue, keep it
        elif residue.id[0] == " ":
            if residue.get_resname() in CANONICAL_RESIDUES:
                return 1
        else:
            return 0


def download_and_uncompress_files(folder, csv_file):
    # Create the input_folder if it doesn't exist
    os.makedirs(folder, exist_ok=True)

    with open(csv_file) as f:
        for line in f:
            id = line.split(',')[0]
            short_id = id[:4]
            file = f"pdb{short_id}.ent.gz"
            # Download the sequence directly from PDB as a compressed file
            subprocess.run(["wget", "-P", folder, f"ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/{file}"])  # noqa
            # The file is uncompressed with gunzip
            subprocess.run(["gunzip", os.path.join(folder, file)])


def rename_files(directory):
    # This function changes all files with a .ent termination to .pdb
    for file in glob.glob(os.path.join(directory, "*.ent")):
        filename = os.path.basename(file)
        new_name = filename[3:7] + ".pdb"
        # Files are saved before cleanup in the raw_PDB folder
        os.rename(file, os.path.join(directory, new_name))


def clean_all_pdb_files(folder):
    # Create the output_folder if it doesn't exist
    os.makedirs(folder, exist_ok=True)

    for file_name in os.listdir(folder):
        if file_name.endswith(".pdb"):
            input_file = os.path.join(folder, file_name)
            # This parses through the .pdb file
            parser = PDBParser()
            structure = parser.get_structure('PDB', input_file)
            io = PDBIO()
            io.set_structure(structure)
            output_file = os.path.join(folder, file_name)
            # Save the NotHetero or main part of the protein as an output file
            io.save(output_file, NotHetero())


def main():
    parser = argparse.ArgumentParser(
                description='Use with CSV file.',
                prog='print_fires.py')

    parser.add_argument('--output_folder',
                        type=str,
                        help='Directory to save the cleaned files',
                        required=True)

    parser.add_argument('--csv_file',
                        type=str,
                        help='Add your csv file',
                        required=True)

    args = parser.parse_args()

    folder = args.output_folder
    csv_file = args.csv_file

    download_and_uncompress_files(folder, csv_file)
    rename_files(folder)
    clean_all_pdb_files(folder)


if __name__ == "__main__":
    main()
