"""
This script is used for downloading, renaming,
and cleaning protein data bank (PDB) files.

Usage:
    python print_fires.py --output_folder <output_folder> --csv_file <csv_file>

Arguments:
    --output_folder (str): Directory to save the cleaned files (required)
    --csv_file (str): Path to the CSV file containing the list of files to download (required)

Example:
    python print_fires.py --output_folder output --csv_file data.csv
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


def download_files(folder, csv_file):
    # Create the input_folder if it doesn't exist
    os.makedirs(folder, exist_ok=True)

    with open(csv_file) as f:
        for line in f:
            id = line.split(',')[0]
            short_id = id[:4]
            file = f"pdb{short_id}.ent.gz"
            # Download the sequence directly from PDB as a compressed file
            subprocess.run(["wget", "-P", folder, f"ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/{file}"])  # noqa


def uncompress_files(directory):
    # This function uncompresses all .gz files in the directory
    for file in glob.glob(os.path.join(directory, "*.gz")):
        subprocess.run(["gunzip", file])


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
    try:
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

        try:
            download_files(folder, csv_file)
        except Exception as e:
            print(f"An error occurred while downloading files: {str(e)}")

        try:
            uncompress_files(folder)
        except Exception as e:
            print(f"An error occurred while uncompressing files: {str(e)}")

        try:
            rename_files(folder)
        except Exception as e:
            print(f"An error occurred while renaming files: {str(e)}")

        try:
            clean_all_pdb_files(folder)
        except Exception as e:
            print(f"An error occurred while cleaning PDB files: {str(e)}")

    except Exception as e:
        print(f"An error occurred in the main function: {str(e)}")


if __name__ == "__main__":
    main()