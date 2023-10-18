import os
from Bio.PDB import PDBParser, PDBList


def find_max_index(residue_insertion_list, end_residue):
    for index, residue in enumerate(reversed(residue_insertion_list)):
        if residue == end_residue:
            return index


def get_d3to1(aa_list):
    seq = []
    d3to1= {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    for aa in aa_list:
        seq.append(d3to1[aa])
    return ''.join(seq)


def get_residue_list_fromDF(atom_df):
    residue_str_list = [str(x) for x in atom_df['residue_number']]
    residue_insertion_list = [a + b for a, b in zip(residue_str_list, atom_df['insertion'])]
    return residue_insertion_list


def get_residue_set_fromList(residue_insertion_list):
    residue_insertion_list_set = []
    [residue_insertion_list_set.append(x) for x in residue_insertion_list if x not in residue_insertion_list_set]
    return residue_insertion_list_set


def get_resolution(pdb_id):
    pdb_list = PDBList()
    parser = PDBParser(QUIET=True)
    try:
        pdb_list.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')
        structure = parser.get_structure(pdb_id, "pdb"+pdb_id+".ent")
    except FileNotFoundError:
        return None
    resolution = structure.header.get('resolution')
    os.remove("pdb"+pdb_id+".ent")
    return resolution


def checkDir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    return dir_path


def dropDup(x):
  return list(dict.fromkeys(x))


def get_list_contains_str(lst, string):
    return [val for val in lst if string in val]
