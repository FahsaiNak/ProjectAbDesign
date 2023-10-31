import sys
from Bio.PDB import PDBParser
sys.path.insert(0, 'src')
import processing_utils as ut   # noqa
import structure_utils as stut   # noqa


def get_pdb_name(matchfile):
    """
    Extract the PDB name from the match file path.

    Args:
        matchfile (str): The path to the match file.

    Returns:
        str: The extracted PDB name.
    """
    info = ut.get_lineinlist(matchfile)
    return info[0].split(".pdb")[0].split("/")[-1]


def get_all_residues(pdbfile):
    """
    Get a list of all residues in a PDB structure.

    Args:
        pdbfile (str): Path to the PDB file.

    Returns:
        list: A list of residue information (e.g., '5|A|SER')
              for all residues in the structure.
        None if the file is not found.
    """
    residue_list = []
    structure = stut.get_structurefromfile(pdbfile)
    if structure is None:
        return None
    for model in structure:
        for chain in model:
            for residue in chain:
                res = f"{residue.id[1]}|{chain.id}|{residue.resname}"
                residue_list.append(res)
    return residue_list
