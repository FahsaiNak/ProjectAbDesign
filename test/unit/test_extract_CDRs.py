import unittest
import random
import sys
sys.path.insert(0, '../../src')  # noqa
import extract_CDRs_from_PDB
import pandas as pd
from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser


class TestExtractCDR(unittest.TestCase):

    def setUp(self):
        self.chothia_pdb_path = '../Datasets/chothia_pdb_files/'
        self.CDR_pdb_path = '../Datasets/CDR_pdb_files/'
        self.pdb_id = '1f58'
        self.chothia_full_path = self.chothia_pdb_path + self.pdb_id + '.pdb'
        self.pdb_df = PandasPdb().read_pdb(self.chothia_full_path) 

        
    def test_get_letter_to_extract(self):
        
        chain = 'H'
        letter_to_extract = extract_CDRs_from_PDB.get_letter_to_extract(self.pdb_df.df['OTHERS'], chain)
        self.assertEqual(letter_to_extract, 'H')

        chain = 'P'
        letter_to_extract = extract_CDRs_from_PDB.get_letter_to_extract(self.pdb_df.df['OTHERS'], chain)
        self.assertEqual(letter_to_extract, None)

    def test_extract_CDR(self):
    
        chain_pdb_letter = 'H'
        atom_df = self.pdb_df.df['ATOM']
        CDR_bounds = [3,7]
        chain_df = atom_df[atom_df.chain_id == chain_pdb_letter].reset_index(drop=True)

        actual_length = 151-111+1
        CDR = extract_CDRs_from_PDB.extract_CDR(chain_df, CDR_bounds)
        self.assertEqual(len(CDR),actual_length)


        CDR_bounds = [3,900]
        CDR = extract_CDRs_from_PDB.extract_CDR(chain_df, CDR_bounds)
        self.assertEqual(CDR,None)


        
    #     self.assertEqual(get_synchro_data.get_data(self.path_to_data + '/' +self.empty_file_name), [])
    #     self.assertRaises(FileNotFoundError, get_synchro_data.get_data,
    #                       self.path_to_data + '/' +self.dne_file_name)
    
    # def test_get_query_lines(self):
        
    #     country = 'Albania'
    #     num_lines = len(get_synchro_data.get_data(self.path_to_data + '/' + self.fire_file_name,
    #                                                query_value=country,
    #                                                query_col=0))
    #     num_lines = len(get_synchro_data.get_data(self.path_to_data + '/' + self.gdp_file_name))
    #     self.assertEqual(199, num_lines)
                        
    #     num_lines = len(get_synchro_data.get_data(self.path_to_data + '/' + self.empty_file_name))
    #     self.assertEqual(0, num_lines)
        

if __name__ == '__main__':
    unittest.main()

        