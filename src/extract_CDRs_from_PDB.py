import pandas as pd
from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser
import processing_utils
import os
import yaml

# GLOBAL VARIABLES

# set paths
with open('../run/config.yaml', 'r') as yaml_file:
    config_data = yaml.safe_load(yaml_file)
CDR_FILE_PATH = config_data["PDBCDRs"]
CHOTHIA_PDB_FILE_PATH = config_data["PDB_Ab_Chothia"]

#  3 and 1 letter amino acid codes
d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

#  Chothia numbering
CDR_Hpos = {'H1': [26, 32], 'H2': [52, 56], 'H3': [95, 102]}
CDR_Lpos = {'L1': [24, 34], 'L2': [50, 56], 'L3': [89, 97]}

# Chain Letters for heavy and light chain
CHAIN_LETTERS = ['H', 'L']

# max CDR-sequence length = 20
# set to minimize compute ime
MAX_CDR_LENGTH = 20

#  check to make sure Nitrogen is in atom list in pdb files
#  otherwise the files are incomplete
ELEMENT_TO_CHECK = 'N'

# pdb id char num
PDB_ID_CHAR = 4


def get_letter_to_extract(head_df=None, query_value=None):
    '''Given the head data frame from a pdb file,
       return the letter corresponding
       to the chain of interest
        ----------
        head_df: DataFrame
                  a data frame information on what
                  letters are used to represent
                  chans in the pdb file atom_df
        query_value: str
                     the letter corresponding to the
                     chain we would like to extract
        Returns
        -------
        let_extract: str
                           the letter by which the heavy
                           or light chain is labeled
    '''

    head_df.drop_duplicates(subset=['entry'], inplace=True)
    chain_indices = head_df["entry"].str.find(query_value+'CHAIN=') + \
        len(query_value+'CHAIN=')
    first_index = (chain_indices > len(query_value+'CHAIN=')).idxmax()
    if first_index == 0:
        return None
    first_value = chain_indices[first_index]
    let_extract = head_df['entry'].iloc[first_index][first_value:first_value+1]

    return let_extract


def extract_CDR(chain_df=None, CDR_bounds=None):
    '''Extracts atoms, residues, and positions from data frame
        containing the chain of interest for specific CDR residue numbers
        according to Chothia numbering
        Parameters
        ----------
        chain_df: DataFrame
                  a data frame containing atom, residue, chain, and positional
                  information
        CDR_bounds: list
                    first index contains starting point for extraction
                    second index contains ending point for extraction
                    e.g., The CDR L1 bound are [24,34]
        Returns
        -------
        CDR_df: DataFrame
                a data frame containing only the rows in chain_df specified
                by the CDR_bounds
    '''

    start_residue = CDR_bounds[0]
    end_residue = CDR_bounds[1]

    start_residue_df = chain_df[chain_df.residue_number == start_residue]
    if len(start_residue_df) == 0:
        print(f'Incomplete CDR chain in file')
        return None
    start_residue_atom_idx = start_residue_df.index.values[0]
    end_residue_df = chain_df[chain_df.residue_number == end_residue]
    if len(end_residue_df) == 0:
        print(f'Incomplete CDR chain in file')
        return None
    end_residue_atom_idx = end_residue_df.index.values[-1]

    CDR_df = chain_df.iloc[start_residue_atom_idx:end_residue_atom_idx+1]
    CDR_df = CDR_df.reset_index(drop=True)
    return CDR_df


def get_CDR_length(pdb_id=None, CDR_save_name=None):
    '''Returns the number of residues contained within a pdb file of a CDR
        ----------
        pdb_id: str
                corresponds to the pdb file we will load in
        CDR_save_name: str
                       the path and file for the pdb entry to read in
        Returns
        -------
        CDR_length: int
                    the length of the CDR
    '''
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, CDR_save_name)

        for model in structure:

            for chain in model:
                seq = []
                for residue in chain:
                    seq.append(d3to1[residue.resname])
                CDR_seq = ''.join(seq)
        CDR_length = len(CDR_seq)
        return CDR_length
    except Exception:
        print(f'PDB Const Error for file: {CDR_save_name}')
        return None


def save_CDR(CDR_df=None, pdb_id=None, CDR_type=None):
    '''Merges newly extracted atomistic CDR data frame with original pdb file
       and saves this new pdb entry to a specified path
        ----------
        CDR_df: DataFrame
                the atomistic CDR data frame to save
        pdb_id: str
                corresponds to the pdb file from which the CDR_df was extracted
        CDR_type: str
                  the CDR type that this CDR data frame is
    '''

    CDR_save_name = CDR_FILE_PATH + '/' + pdb_id + '_' + CDR_type + '.pdb'

    pdb_to_save = PandasPdb().read_pdb(CHOTHIA_PDB_FILE_PATH+'/'+pdb_id+'.pdb')
    pdb_to_save.df['ATOM'] = pd.merge(pdb_to_save.df['ATOM'],
                                      CDR_df,
                                      how='inner')
    pdb_to_save.to_pdb(path=CDR_save_name,
                       records=['ATOM', 'OTHERS', 'ANISOU'],
                       gz=False,
                       append_newline=True)


def main():

    processing_utils.checkDir(CHOTHIA_PDB_FILE_PATH)
    processing_utils.checkDir(CDR_FILE_PATH)

    files = os.listdir(CHOTHIA_PDB_FILE_PATH)

    for index, file in enumerate(files):

        pdb_id = file[:-PDB_ID_CHAR]
        print('Index:', index, pdb_id)
        pdb_df = PandasPdb().read_pdb(CHOTHIA_PDB_FILE_PATH + '/' + file)

        head_df = pdb_df.df['OTHERS']
        atom_df = pdb_df.df['ATOM']

        atom_name_series = atom_df['atom_name']

        # removes pdb files that do not contain full atom list
        # by checking for presence of nitrogen atom
        # assumes that if N is present then it's full chain
        if ELEMENT_TO_CHECK in atom_name_series.values:

            for chain_letter in CHAIN_LETTERS:

                chain_pdb_letter = get_letter_to_extract(
                    head_df=head_df,
                    chain_letter=chain_letter
                )

                if chain_pdb_letter is None:
                    print(f'No H chain in file')
                else:

                    chain_df = atom_df[atom_df.chain_id == chain_pdb_letter]
                    chain_df = chain_df.reset_index(drop=True)
                    for CDR_type, CDR_bounds in CDR_Hpos.items():

                        CDR_df = extract_CDR(chain_df=chain_df,
                                             CDR_bounds=CDR_bounds)

                        save_CDR(CDR_df=CDR_df,
                                 pdb_id=pdb_id,
                                 CDR_type=CDR_type)

                        CDR_save_name = CDR_FILE_PATH + '/' + pdb_id \
                            + '_' + CDR_type + '.pdb'
                        CDR_length = get_CDR_length(
                            pdb_id=pdb_id,
                            CDR_save_name=CDR_save_name
                        )

                        if CDR_length is None:
                            os.remove(CDR_save_name)
                        elif CDR_length >= MAX_CDR_LENGTH:
                            print(f'{CDR_save_name} too long. Removing file')
                            os.remove(CDR_save_name)
        else:
            print(f"{ELEMENT_TO_CHECK} is not in the DataFrame.")
            os.remove(CHOTHIA_PDB_FILE_PATH + '/' + file)
        print('----------------------------------------------')


if __name__ == '__main__':
    main()
