import unittest
import sys
sys.path.insert(0, '../../src')  # noqa
import fragment_CDRs
from biopandas.pdb import PandasPdb


class TestExtractCDR(unittest.TestCase):

    def test_slide_window(self):

        frag_dict = fragment_CDRs.slide_window('1f58_H1.pdb')
        self.assertEqual(frag_dict['sequence'][0], 'GYSI')


if __name__ == '__main__':
    unittest.main()
