import unittest
import sys
sys.path.insert(0, '../../src')  # noqa
import extract_CDRs_from_PDB
from biopandas.pdb import PandasPdb


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
        self.assertIsNone(letter_to_extract)

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
        self.assertIsNone(CDR)

    
    def test_get_CDR_length(self):
        
        length = extract_CDRs_from_PDB.get_CDR_length(self.pdb_id,
                                                      self.CDR_pdb_path + self.pdb_id + '_H1.pdb')

        actual_length = 32-26 + 1 + 4
        self.assertEqual(length,actual_length)

        length = extract_CDRs_from_PDB.get_CDR_length('bubba',
                                                      self.CDR_pdb_path + self.pdb_id + '_H1.pdb')

        self.assertEqual(length,actual_length)

        length = extract_CDRs_from_PDB.get_CDR_length('bubba',
                                                      self.CDR_pdb_path + self.pdb_id + '.pdb')

        self.assertIsNone(length)        

if __name__ == '__main__':
    unittest.main()

        