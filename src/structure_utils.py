import os
from Bio.PDB import PDBParser, PDBList
import master_utils as mut
import processing_utils as ut


def find_max_index(residue_insertion_list, end_residue):
    '''Takes in a list of residues w/ insertions and returns
        the max index for an end residue
        Parmams
        -------
        residue_insertion_list: Series
                                a series of residues w/ insertions
        end_residue: str
                     the end residue to find the index of
                     e.g. '30' or '30A'
        Returns
        ------
        index: int
               the maximum index in the series with the end_residue
        '''
    for index, residue in enumerate(reversed(residue_insertion_list)):
        if residue == end_residue:
            return index


def get_residue_list_fromDF(atom_df):
    residue_str_list = [str(x) for x in atom_df['residue_number']]
    residue_insertion_list = [a + b for a, b in zip(residue_str_list, atom_df['insertion'])]
    return residue_insertion_list


def get_residue_set_fromList(residue_insertion_list):
    residue_insertion_list_set = []
    [residue_insertion_list_set.append(x) for x in residue_insertion_list if x not in residue_insertion_list_set]
    return residue_insertion_list_set


def get_structurefrompdb(pdb_id):
    pdb_list = PDBList()
    parser = PDBParser(QUIET=True)
    try:
        pdb_list.download_pdb_files([pdb_id], pdir='.', file_format='pdb')
        structure = parser.get_structure(pdb_id, "pdb"+pdb_id+".ent")
    except FileNotFoundError:
        return None
    os.remove("pdb"+pdb_id+".ent")
    return structure


def get_structurefromfile(pdbfile):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(mut.get_pdb_name(pdbfile), pdbfile)
    except FileNotFoundError:
        return None
    return structure


def get_resolution(pdb_id):
    structure = get_structurefrompdb(pdb_id)
    resolution = structure.header.get('resolution')
    return resolution


def find_SASAcontact(full_struct, crop_struc):
    agsasa_list = []
    for chain in crop_struc[0]:
        for res in chain:
            try:
                sasa_val_crop = round(res.sasa, 2)
                sasa_val_full = round(full_struct[0][chain.id][res.id].sasa, 2)
            except AttributeError:
                continue
            if sasa_val_crop >= sasa_val_full+0.4:
                if res.resname in ut.call_d3to1().keys():
                    res_out = f"{res.id[1]}|{chain.id}|{res.resname}"
                    agsasa_list.append(res_out)
    return agsasa_list
