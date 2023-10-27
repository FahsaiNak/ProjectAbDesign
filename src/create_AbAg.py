import argparse
import pandas as pd
import os
import pickle
from glob import glob

def get_args():
    parser = argparse.ArgumentParser(
        description='combine Ab-like, Ag-like information for all',
        prog='create_AbAg')
    parser.add_argument('--dir', type=str,
                        help='directory that stores each pickle file of Ab-like and Ag-like',
                        required=True)
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    datadir = args.dir
    output_df = pd.concat([pd.read_pickle(os.path.join(datadir, x)) for x in os.listdir(datadir)])
    output_df.reset_index(drop=True, inplace=True)
    output_df.to_pickle(os.path.join('/'.join(datadir.split("/")[:-1]), "AbAg.pkl"))
    for filename in glob(os.path.join(datadir, "*.pkl")):
        os.remove(filename)
    os.rmdir(datadir)


if __name__ == '__main__':
    main()
