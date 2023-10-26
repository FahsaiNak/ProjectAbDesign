import sys
import argparse
import pandas as pd
from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser, PDBList
import processing_utils as ut
import structure_utils as frut
import os


def get_args():
    parser = argparse.ArgumentParser(
        description='Sub-fragmenting CDRs by a sliding window method',
        prog='CDR_fragmentation')
    parser.add_argument('--path', type=str,
                        help='Path to retrieve CDR dataset',
                        required=True)
    parser.add_argument('--savePath', type=str,
                        help='Path to save fragmented-CDRs',
                        required=True)
    args = parser.parse_args()
    return args


def getInfo():
    args = get_args()
    files = os.listdir(args.path)
    pdb_lst = ut.dropDup([_[:4] for _ in files])
    result_dict = {"name":[], "resolution":[], "sequence":[],
                   "start_residue_idx":[], "max_index_from_end":[]}
    for i, pdb in enumerate(pdb_lst):
        print(i+1, "/", len(pdb_lst), pdb, "is processing")
        resolution = frut.get_resolution(pdb)
        cdr_files = ut.get_list_contains_str(files, pdb)
        for file in cdr_files:
            try:
                CDR_pdb = PandasPdb().read_pdb(os.path.join(args.path, file))
                atom_df = CDR_pdb.df["ATOM"]
                residue_insertion_list = frut.get_residue_list_fromDF(atom_df)
                residue_insertion_list_set = frut.get_residue_set_fromList(residue_insertion_list)
                num_residues = len(residue_insertion_list_set)
                for start_residue_number in range(0,len(residue_insertion_list_set)-4+1):
                    start_residue = residue_insertion_list_set[start_residue_number]
                    start_residue_idx = residue_insertion_list.index(start_residue)
                    for fragment_length in range(4, num_residues-start_residue_number+1):
                        end_residue = residue_insertion_list_set[start_residue_number+fragment_length-1]
                        max_index_from_end = frut.find_max_index(residue_insertion_list, end_residue)
                        if max_index_from_end == 0:
                            max_index_from_end = -len(residue_insertion_list)
                        fragment_df = atom_df.iloc[start_residue_idx:-max_index_from_end]
                        aa_list = fragment_df.drop_duplicates(subset=["residue_number", "insertion"]).residue_name.to_list()
                        seq = ut.get_d3to1(aa_list)
                        result_dict["name"].append(file)
                        result_dict["resolution"].append(resolution)
                        result_dict["sequence"].append(seq)
                        result_dict["start_residue_idx"].append(start_residue_idx)
                        result_dict["max_index_from_end"].append(max_index_from_end)
            except ValueError:
                print(f'{file} not pdb')
                continue
        print(len(cdr_files), "CDRs have been processed")
    return result_dict


def main():
    args = get_args()
    result_dict = getInfo()
    result_df = pd.DataFrame(result_dict)
    result_df.sort_values(by=["sequence", "resolution"], inplace=True)
    result_df.drop_duplicates(subset=["sequence"], keep='first', inplace=True)
    result_df.reset_index(drop=True, inplace=True)
    for ind in result_df.index:
        file = result_df.loc[ind, "name"]
        start_residue_idx = result_df.loc[ind, "start_residue_idx"]
        max_index_from_end = result_df.loc[ind, "max_index_from_end"]
        pdb_to_save = PandasPdb().read_pdb(os.path.join(args.path, file))
        atom_df = pdb_to_save.df["ATOM"]
        fragment_df = atom_df.iloc[start_residue_idx:-max_index_from_end]
        residue_list = frut.get_residue_list_fromDF(fragment_df)
        pdb_to_save.df['ATOM'] = pd.merge(pdb_to_save.df['ATOM'], fragment_df, how='inner')
        pdb_save_file_name = os.path.join(ut.checkDir(args.savePath), file[:-4]+'_frag_'+residue_list[0]+'_'+residue_list[-1]+'.pdb')
        pdb_to_save.to_pdb(path=pdb_save_file_name, records=['ATOM', 'OTHERS', 'ANISOU'], gz=False, append_newline=True)
        print(ind+1,"/",len(result_df.index), pdb_save_file_name, "saved")


if __name__ == '__main__':
    main()

