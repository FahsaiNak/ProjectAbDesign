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
    parser.add_argument('--indices_only', type=bool,
                        help='If true, will output AbAg regions as indices,'
                             'only, stripping out chain, amino acid code, and '
                             'separators', required=False, default=True)
    args = parser.parse_args()
    return args


def main():
    # Converts AbAg to dictionary with pdb_id as the key, and the Ab and Ag-
    # like regions as the values.
    # Values are formatted as a list of tuples of the form (Ablike, Aglike)
    args = get_args()
    abag = ut.open_pickle(args.abag_filename)
    ids_list = list(abag['pdb_id'])
    ab_list = list(abag['ablike'])
    ag_list = list(abag['aglike'])
    # Extracting indices from regions and exporting as .csv
    # Reformat strings
    ab_list_parsed = [ut.parse_string_to_list(x) for x in ab_list]
    ag_list_parsed = [ut.parse_string_to_list(x) for x in ag_list]

    # Below is the dictionary with unabridged regions
    abag_dict_full_strings = ut.refactor_abag(ids_list, ab_list, ag_list)
    # Below is the dictionary with regions parsed to only include AA indices
    abag_dict_indices_only = ut.refactor_abag(ids_list, ab_list_parsed,
                                              ag_list_parsed)
    # I don't know what to do with these dictionaries! the export format that
    # makes the most sense to me/that we'll use for the visualization is below.
    # This format is one row per ag-ab pair, with some pdb_ids having multiple
    # rows
    if args.indices_only:
        export_df = pd.DataFrame(zip(ids_list, ab_list_parsed, ag_list_parsed))
    else:
        export_df = pd.DataFrame(zip(ids_list, ab_list, ag_list))
    export_df.columns = ['pdb_id', 'ab_like', 'ag_like']
    export_df.to_csv(args.output_filename, index=False)


if __name__ == '__main__':
    main()