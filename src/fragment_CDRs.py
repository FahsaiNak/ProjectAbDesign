import sys
import argparse
import pandas as pd
from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser, PDBList
import processing_utils as ut
import structure_utils as frut #TODO change from frut
import os
import yaml

#GLOBAL VARIABLES

# Set paths for file
with open('../run/config.yaml', 'r') as yaml_file:
    config_data = yaml.safe_load(yaml_file)
CDR_FILE_PATH = config_data["PDBCDRs"]
CDR_FRAGMENT_PATH = config_data["PDBCDRfrag"]

# set the minimum fragment length according to Rangel et al.
MIN_FRAG_LENGTH = 4     

def get_CDR_frag_dict():
    '''Performs sliding window on all files in the CDR directory
    Returns
    ------
    fragment_dict: dict
                   a dictionary containing 
                   - CDR_file: the CDR_file used for fragment generation
                   - resolution: the resolution for this CDR_file
                   - sequence: the single letter fragment sequence
                   - start_residue_idx: the beginnning atom index to begin
                     extraction from the CDR atomistic data frame
                   - last_residue_index_from_end: the last atom index from
                     the end for extraction from the CDR atomistic data frame
    '''

    fragment_dict = {"CDR_file":[], "resolution":[], "sequence":[],
                    "start_residue_idx":[], "last_residue_index_from_end":[]}

    CDR_files = os.listdir(CDR_FILE_PATH)
    print(CDR_files)
    pdb_id_lst = ut.dropDup([_[:4] for _ in CDR_files]) 
    
    #TODO I need to figure out how to not make four for loops
    for idx, pdb_id in enumerate(pdb_id_lst):
        
        print(idx+1, "/", len(pdb_id_lst), pdb_id, "is processing")
        resolution = frut.get_resolution(pdb_id)
        CDR_files = ut.get_list_contains_str(CDR_files, pdb_id)
        for CDR_file in CDR_files:

            next_frag_dict = slide_window(CDR_file=CDR_file, resolution=resolution)

            fragment_dict["CDR_file"].extend(next_frag_dict["CDR_file"])
            fragment_dict["resolution"].extend(next_frag_dict["resolution"])
            fragment_dict["sequence"].extend(next_frag_dict["sequence"])
            fragment_dict["start_residue_idx"].extend(next_frag_dict["start_residue_idx"])
            fragment_dict["last_residue_index_from_end"].extend(next_frag_dict["last_residue_index_from_end"])

    return fragment_dict

def slide_window(CDR_file = None, resolution=None):

    fragment_dict = {"CDR_file":[], "resolution":[], "sequence":[],
                    "start_residue_idx":[], "last_residue_index_from_end":[]}

    try:

        CDR_pdb = PandasPdb().read_pdb(os.path.join(CDR_FILE_PATH, CDR_file))
        atom_df = CDR_pdb.df["ATOM"]
        residue_insertion_list = frut.get_residue_list_fromDF(atom_df)
        residue_insertion_list_set = frut.get_residue_set_fromList(residue_insertion_list)
        num_residues = len(residue_insertion_list_set)

        #perform sliding window on atom data frame
        for start_residue_number in range(0, len(residue_insertion_list_set) - MIN_FRAG_LENGTH + 1):

            start_residue = residue_insertion_list_set[start_residue_number]
            start_residue_idx = residue_insertion_list.index(start_residue)

            for fragment_length in range(MIN_FRAG_LENGTH, num_residues - start_residue_number + 1):

                last_residue = residue_insertion_list_set[start_residue_number+fragment_length - 1]
                last_residue_index_from_end = frut.find_max_index(residue_insertion_list, last_residue)

                if last_residue_index_from_end == 0:
                    last_residue_index_from_end = -len(residue_insertion_list)

                fragment_df = atom_df.iloc[start_residue_idx:-last_residue_index_from_end]

                aa_list = fragment_df.drop_duplicates(subset=["residue_number", "insertion"]).residue_name.to_list()
                seq = ut.get_d3to1(aa_list)
                fragment_dict["CDR_file"].append(CDR_file)
                fragment_dict["resolution"].append(resolution)
                fragment_dict["sequence"].append(seq)
                fragment_dict["start_residue_idx"].append(start_residue_idx)
                fragment_dict["last_residue_index_from_end"].append(last_residue_index_from_end)

    except ValueError:

        print(f'{CDR_file} not pdb')

    return fragment_dict

def main():

    fragment_dict = get_CDR_frag_dict()
    result_df = pd.DataFrame(fragment_dict)
    result_df.sort_values(by=["sequence", "resolution"], inplace=True)
    result_df.drop_duplicates(subset=["sequence"], keep='first', inplace=True)
    result_df.reset_index(drop=True, inplace=True)


    for ind in result_df.index:
        file = result_df.loc[ind, "CDR_file"]
        start_residue_idx = result_df.loc[ind, "start_residue_idx"]
        max_index_from_end = result_df.loc[ind, "last_residue_index_from_end"]
        
        pdb_to_save = PandasPdb().read_pdb(os.path.join(CDR_FILE_PATH, file))
        atom_df = pdb_to_save.df["ATOM"]
        fragment_df = atom_df.iloc[start_residue_idx:-max_index_from_end]
        residue_list = frut.get_residue_list_fromDF(fragment_df)
        pdb_to_save.df['ATOM'] = pd.merge(pdb_to_save.df['ATOM'], fragment_df, how='inner')
        pdb_save_file_name = os.path.join(ut.checkDir(CDR_FRAGMENT_PATH),
                                         file[:-4]+'_frag_'+residue_list[0]+'_'+residue_list[-1]+'.pdb')
        pdb_to_save.to_pdb(path=pdb_save_file_name, records=['ATOM', 'OTHERS', 'ANISOU'],
                           gz=False, append_newline=True)
        print(ind+1,"/",len(result_df.index), pdb_save_file_name, "saved")


if __name__ == '__main__':
    main()

