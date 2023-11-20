import unittest
import sys
sys.path.insert(0, '../../src')  # noqa
import fragment_CDRs
from biopandas.pdb import PandasPdb


class TestExtractCDR(unittest.TestCase):

    def test_get_CDR_frag_dict(self):

        fragment_dict = fragment_CDRs.get_CDR_frag_dict()
        actual_frag_seq = 'SQGV'
        self.assertEqual(actual_frag_seq, fragment_dict['sequence'][0])

    def test_slide_window(self):

        frag_dict = fragment_CDRs.slide_window('1f58_H1.pdb')
        self.assertEqual(frag_dict['sequence'][0], 'SQGV')


if __name__ == '__main__':
    unittest.main()
