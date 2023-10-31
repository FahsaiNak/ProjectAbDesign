import sys
import argparse
from glob import glob
import os
import pickle
import pandas as pd
sys.path.insert(0, 'src')
import master_utils as mut  # noqa


def get_args():
    parser = argparse.ArgumentParser(
        description='Get CDR-like information in PDB structure',
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
    """
    Main function to process matched PDB structures
    and extract CDR-like information.
    """
    args = get_args()
    pdbid = args.pdb
    datadir = args.data_dir

    # Initialize a dictionary to store CDR-like information.
    collection = {"ab_pos1": [], "ab_pos2": [], "ab_all": []}

    # Check if there are any matching structure files in the data directory.
    if len(glob(os.path.join(datadir, "*."+pdbid+".struct"))) == 0:
        sys.exit(1)

    # Loop through directories containing matching structures.
    for dir_path in glob(os.path.join(datadir, "*."+pdbid+".struct")):
        print(dir_path)

        # Loop through PDB files within each directory.
        for match_path in glob(os.path.join(dir_path, "*.pdb")):
            res_all = mut.get_all_residues(match_path)
            out = ','.join(res_all)
            info_lst = [res_all[0], res_all[-1], out]

            # Populate the collection dictionary.
            for k, v in zip(collection.keys(), info_lst):
                collection[k].append(v)

    # Create a Pandas DataFrame from the collected information.
    df = pd.DataFrame(collection)

    # Sort and remove duplicates from the DataFrame.
    df.sort_values(by=["ab_pos1", "ab_pos2", "ab_all"],
                   ascending=[True, False, False], inplace=True)
    df.drop_duplicates(subset=["ab_pos1", "ab_pos2", "ab_all"], inplace=True)
    df.drop_duplicates(subset=["ab_pos1"], inplace=True)
    df.drop_duplicates(subset=["ab_pos2"], inplace=True)
    df.reset_index(drop=True, inplace=True)

    # Save the "ab_all" column as a Pickle file.
    df["ab_all"].to_pickle(os.path.join(datadir, pdbid+".pkl"))


if __name__ == '__main__':
    main()
