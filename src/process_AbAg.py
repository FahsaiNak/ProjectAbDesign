import pandas as pd
import pickle
import argparse
import process_AbAg_utils as ut


def get_args():
    parser = argparse.ArgumentParser(
        description='Get CDR-like information in PDB structure',
        prog='process_AbAg')
    parser.add_argument('--abag_filename', type=str,
                        help='filename of abag file, in pickle format',
                        required=True)
    parser.add_argument('--output_filename', type=str,
                        help='name of output file',required=True)
    args = parser.parse_args()
    return args


def main():
    # Converts AbAg in pickle format to csv with pdb id in the first column,
    # ab-like region indices in the second column, ab-like chain letters in the
    # third column, ag-like region indices in the fourth column, and ag-like
    # chain letters in the fifth column
    args = get_args()
    abag = ut.open_pickle(args.abag_filename)
    ids_list = list(abag['pdb_id'])
    ab_list = list(abag['ablike'])
    ag_list = list(abag['aglike'])
    # Extracting indices from regions and exporting as .csv
    # Reformat strings
    ab_list_ints = [ut.parse_string_to_list_ints(x) for x in ab_list]
    ag_list_ints = [ut.parse_string_to_list_ints(x) for x in ag_list]
    ab_list_chains = [ut.parse_string_to_list_chain_letters(x) for x in
                      ab_list]
    ag_list_chains = [ut.parse_string_to_list_chain_letters(x) for x in
                      ag_list]
    export_df = pd.DataFrame(zip(ids_list, ab_list_ints, ab_list_chains,
                                 ag_list_ints, ag_list_chains))
    export_df.columns = ['pdb_id', 'ab-like integers', 'ab-like chains',
                         'ag-like integers', 'ag-like chains']
    export_df.to_csv(args.output_filename, index=False)


if __name__ == '__main__':
    main()