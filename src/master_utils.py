import processing_utils as ut
from Bio.PDB import PDBParser


def get_pdb_name(matchfile):
    info = ut.get_lineinlist(matchfile)
    return info[0].split(".pdb")[0].split("/")[-1]


def get_residue_no(matchfile, position):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(get_pdb_name(matchfile), matchfile)
    residue_list = [_ for _ in structure.get_residues()]
    res_no = residue_list[position].id[1]
    res_chain = residue_list[position].get_parent().id
    return [res_chain, res_no]


def get_seq(matchfile):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(get_pdb_name(matchfile), matchfile)
    resname_list = [_.resname for _ in structure.get_residues()]
    return ut.get_d3to1(resname_list)


