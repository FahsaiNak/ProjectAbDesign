import sys
import argparse
import pandas as pd
import os
import pickle
from glob import glob


def get_args():
    parser = argparse.ArgumentParser(
        description='Combine Ab-like, Ag-like information for all',
        prog='create_AbAg')
    parser.add_argument('--dir', type=str,
                        help='Directory that stores each\
                        pickle file of Ab-like and Ag-like',
                        required=True)
    args = parser.parse_args()
    return args


def main():
    """
    Main function to combine Ab-like and Ag-like information
    from multiple files into a single file.
    """
    args = get_args()
    datadir = args.dir

    # Get a list of all AbAg.pkl files in the specified directory.
    files = glob(os.path.join(datadir, "*.AbAg.pkl"))

    if len(files) == 0:
        sys.exit(1)

    # Read and concatenate data from all AbAg.pkl files.
    output_df = pd.concat([pd.read_pickle(file) for file in files])

    # Reset the index in the concatenated DataFrame.
    output_df.reset_index(drop=True, inplace=True)

    # Save the combined data as AbAg.pkl in the parent directory
    # of the specified directory.
    output_df.to_pickle(os.path.join('/'.join(datadir.split("/")[:-1]),
                                     "AbAg.pkl"))


if __name__ == '__main__':
    main()
