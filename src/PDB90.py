import os
import time
import glob
import subprocess
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
    # by removing small molecules or ions that are not a main part of the protein
    def accept_residue(self, residue):
        # If the sequence contains water molecules (HOH or WAT), return 0, removing them
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

def download_and_uncompress_files(folder, csv_directory):
    # Create the input_folder if it doesn't exist
    os.makedirs(folder, exist_ok=True)
    
    # somePDB.csv is a txt file that contains a 4 alphanumerical ID for each sequence
    # somePDB only contains the IDs for 10 sequences, the representativePDB.csv containd the IDs for 57148 sequences
    # Running the script with representativePDB.csv takes several hours or even days
    # Run only with that file if you need to construct the full database
    with open(os.path.join(csv_directory, 'somePDB.csv')) as f:
        for line in f:
            id = line.split(',')[0]
            short_id = id[:4]
            file = f"pdb{short_id}.ent.gz"
            # This downloads the sequence directly from PDB as a compressed file
            subprocess.run(["wget", "-P", folder, f"ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/{file}"])
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
            # This saves the NotHetero or main part of the protein as an output file
            io.save(output_file, NotHetero())

def main():
    folder = './Datasets/all_PDB'
    csv_directory = '../' 
    
    download_and_uncompress_files(folder, csv_directory)
    rename_files(folder)
    clean_all_pdb_files(folder)

if __name__ == "__main__":
    main()

