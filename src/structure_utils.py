import os
from Bio.PDB import PDBParser, PDBList
sys.path.insert(0, 'src')
import master_utils as mut   # noqa
import processing_utils as ut   # noqa


def find_max_index(residue_insertion_list, end_residue):
    """
    Find the maximum index of a specific end residue
    in a list of residue insertions.

    Args:
        residue_insertion_list (list): List of residue insertions.
        end_residue (str): The end residue to search for.

    Returns:
        int: The maximum index of the end residue
        in the list, or None if not found.
    """
    for index, residue in enumerate(reversed(residue_insertion_list)):
        if residue == end_residue:
            return index


def get_residue_list_fromDF(atom_df):
    """
    Get a list of residue insertions from a DataFrame of atoms.

    Args:
        atom_df (pandas.DataFrame): DataFrame containing atom information,
        including residue numbers and insertions.

    Returns:
        list: List of residue insertions.
    """
    residue_str_list = [str(x) for x in atom_df['residue_number']]
    residue_insertion_list = [a + b for a, b in zip(residue_str_list, atom_df['insertion'])]   # noqa

    return residue_insertion_list


def get_residue_set_fromList(residue_insertion_list):
    """
    Get a set of unique residue insertions from a list of residue insertions.

    Args:
        residue_insertion_list (list): List of residue insertions.

    Returns:
        set: Set of unique residue insertions.
    """
    residue_insertion_list_set = []
    [residue_insertion_list_set.append(x) for x in residue_insertion_list if x not in residue_insertion_list_set]   # noqa
    return residue_insertion_list_set


def get_structurefrompdb(pdb_id):
    """
    Download and parse a PDB structure by its PDB ID.

    Args:
        pdb_id (str): PDB ID of the structure to download.

    Returns:
        Bio.PDB.Structure.Structure: Parsed PDB structure object,
                                     or None if the file is not found.
    """
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
    """
    Parse a PDB structure from a local file.

    Args:
        pdbfile (str): Path to the local PDB file.

    Returns:
        Bio.PDB.Structure.Structure: Parsed PDB structure object,
                                     or None if the file is not found.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(mut.get_pdb_name(pdbfile), pdbfile)
    except FileNotFoundError:
        return None
    return structure


def get_resolution(pdb_id):
    """
    Get the resolution of a PDB structure by its PDB ID.

    Args:
        pdb_id (str): PDB ID of the structure to retrieve resolution data.

    Returns:
        float: The resolution of the structure,
               or None if the file is not found.
    """
    structure = get_structurefrompdb(pdb_id)
    resolution = structure.header.get('resolution')
    return resolution


def find_SASAcontact(full_struct, crop_struc):
    """
    Find SASA (Solvent Accessible Surface Area) contacts
    in the cropped structure that have significantly higher SASA values
    compared to the corresponding residues in the full structure.

    Args:
        full_struct (Bio.PDB.Structure.Structure): The full PDB structure.
        crop_struc (Bio.PDB.Structure.Structure): The cropped PDB structure.

    Returns:
        agsasa_list (list): A list of residues with SASA contacts.
    """
    agsasa_list = []
    for chain in crop_struc[0]:
        for res in chain:
            try:
                sasa_val_crop = round(res.sasa, 2)
                sasa_val_full = round(full_struct[0][chain.id][res.id].sasa, 2)
            except AttributeError:
                continue
            if sasa_val_crop >= sasa_val_full + 0.4:
                if res.resname in ut.call_d3to1().keys():
                    res_out = f"{res.id[1]}|{chain.id}|{res.resname}"
                    agsasa_list.append(res_out)
    return agsasa_list
