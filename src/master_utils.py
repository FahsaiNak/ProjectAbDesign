import processing_utils as ut
import structure_utils as stut
from Bio.PDB import PDBParser


def get_pdb_name(matchfile):
    info = ut.get_lineinlist(matchfile)
    return info[0].split(".pdb")[0].split("/")[-1]


def get_all_residues(pdbfile):
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
