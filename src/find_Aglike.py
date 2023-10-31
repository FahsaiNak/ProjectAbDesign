import sys
import argparse
import pandas as pd
import os
import pickle
from Bio.PDB import ShrakeRupley
sys.path.insert(0, 'src')
import processing_utils as ut  # noqa
import structure_utils as stut  # noqa


def get_args():
    parser = argparse.ArgumentParser(
        description='Find Ag-like structure in PDB structure',
        prog='get_Ag')
    parser.add_argument('--file', type=str,
                        help='Pickle file that stores Ab-like information',
                        required=True)
    args = parser.parse_args()
    return args


def main():
    """
    Main function to find Ag-like structures
    in a PDB file based on Ab-like information.
    """
    result_dict = {"pdb_id": [], "ablike": [], "ablike_seq": [],
                   "aglike": [], "aglike_seq": []}

    sr = ShrakeRupley()
    args = get_args()
    file = args.file
    pdbid = file.split("/")[-1].split(".")[0]
    datadir = "/".join(file.split("/")[:-1])

    # Check if the input pickle file exists.
    if os.path.isfile(file) is False:
        sys.exit(1)

    # Read Ab-like information from the pickle file.
    df = pd.read_pickle(file)

    # Get the full PDB structure for the specified PDB ID.
    full_struct = stut.get_structurefrompdb(pdbid)

    if full_struct is None:
        sys.exit(1)

    # Compute the Shrake-Rupley SASA for the full structure.
    sr.compute(full_struct[0], level="R")
    crop_struc = full_struct.copy()

    for i, abres in enumerate(df):
        sasa_res_list = list()
        abres_list = abres.split(",")
        abres_no = [int(_.split("|")[0]) for _ in abres_list]
        abres_chain = [_.split("|")[1] for _ in abres_list]
        abres_seq = ut.get_d3to1([_.split("|")[2] for _ in abres_list])
        crop_struc = full_struct.copy()

        # Remove Ab residues from the structure.
        for ab_chain, ab_no in zip(abres_chain, abres_no):
            res_id = full_struct[0][ab_chain][ab_no].id
            crop_struc[0][ab_chain].detach_child(res_id)

        # Compute SASA for the modified structure.
        sr.compute(crop_struc[0], level="R")

        # Find SASA contacts that are not part of Ab-like residues.
        sasa_res_out = stut.find_SASAcontact(full_struct, crop_struc)
        sasa_res_list.extend([_ for _ in sasa_res_out if _ not in abres_list])
        agres_list = ut.dropDup(sasa_res_list)
        agres = ",".join(agres_list)
        agres_seq = ut.get_d3to1([_.split("|")[2] for _ in agres_list])

        # Print results.
        print("CDR-like", abres_seq, "Ag-like", agres_seq)

        if len(agres_list) != 0:
            # Store the results in the result_dict.
            result_dict["pdb_id"].append(pdbid)
            result_dict["ablike"].append(abres)
            result_dict["ablike_seq"].append(abres_seq)
            result_dict["aglike"].append(agres)
            result_dict["aglike_seq"].append(agres_seq)

    result_df = pd.DataFrame(result_dict)

    # Remove duplicates and reset the index in the result DataFrame.
    result_df.drop_duplicates(subset=["ablike_seq", "aglike_seq"],
                              inplace=True)
    result_df.reset_index(drop=True, inplace=True)

    # Save the results to a Pickle file.
    result_df[["pdb_id", "ablike", "aglike"]].to_pickle(
        os.path.join(datadir, pdbid+".AbAg.pkl"))


if __name__ == '__main__':
    main()
