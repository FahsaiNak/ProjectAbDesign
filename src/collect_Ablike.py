import sys
import argparse
from glob import glob
import os
import pickle
import pandas as pd
sys.path.insert(0, 'src')
import master_utils as mut


def get_args():
    parser = argparse.ArgumentParser(
        description='get CDR-like information in PDB structure',
        prog='get_Ab')
    parser.add_argument('--pdb', type=str,
                        help='PDB id',
                        required=True)
    parser.add_argument('--data_dir', type=str,
                        help='Directory that stores match structures',
                        required=True)
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    pdbid = args.pdb
    datadir = args.data_dir
    collection = {"ab_pos1":[], "ab_pos2":[], "ab_all":[]}
    if len(glob(os.path.join(datadir, "*."+pdbid+".struct"))) == 0:
        sys.exit(1)
    for dir_path in glob(os.path.join(datadir, "*."+pdbid+".struct")):
        print(dir_path)
        for match_path in glob(os.path.join(dir_path, "*.pdb")):
            res_all = mut.get_all_residues(match_path)
            out = ','.join(res_all)
            info_lst = [res_all[0], res_all[-1], out]
            for k, v in zip(collection.keys(), info_lst):
                collection[k].append(v)
    df = pd.DataFrame(collection)
    df.sort_values(by=["ab_pos1", "ab_pos2", "ab_all"], ascending=[True,False,False], inplace=True)
    df.drop_duplicates(subset=["ab_pos1", "ab_pos2", "ab_all"], inplace=True)
    df.drop_duplicates(subset=["ab_pos1"], inplace=True)
    df.drop_duplicates(subset=["ab_pos2"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    df["ab_all"].to_pickle(os.path.join(datadir,pdbid+".pkl"))


if __name__ == '__main__':
    main()
