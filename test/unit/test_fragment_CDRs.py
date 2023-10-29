import unittest
import sys
sys.path.insert(0, '../../src')  # noqa
import fragment_CDRs
from biopandas.pdb import PandasPdb


class TestExtractCDR(unittest.TestCase):
        
    def test_get_CDR_frag_dict(self):
        
        fragment_dict = fragment_CDRs.get_CDR_frag_dict()
        actual_frag_seq = 'SQGV'
        fragment_dict['sequence'][0]
        self.assertEqual(actual_frag_seq,fragment_dict['sequence'][0])

    
if __name__ == '__main__':
    unittest.main()

        