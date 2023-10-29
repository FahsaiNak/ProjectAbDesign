import unittest
import sys
sys.path.insert(0, '../../src')  # noqa
import fragment_CDRs
from biopandas.pdb import PandasPdb


class TestExtractCDR(unittest.TestCase):

    def setUp(self):
        self.chothia_pdb_path = '../Datasets/chothia_pdb_files/'
        self.CDR_pdb_path = '../Datasets/CDR_pdb_files/'
        self.pdb_id = '1f58'
        self.chothia_full_path = self.chothia_pdb_path + self.pdb_id + '.pdb'
        self.pdb_df = PandasPdb().read_pdb(self.chothia_full_path) 

        
    def test_get_CDR_frag_dict(self):
        
        fragment_dict = fragment_CDRs.get_CDR_frag_dict()
        actual_frag_seq = 'SQGV'
        fragment_dict['sequence'][0]
        self.assertEqual(actual_frag_seq,fragment_dict['sequence'][0])

    
if __name__ == '__main__':
    unittest.main()

        