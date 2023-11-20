import unittest
import sys
sys.path.insert(0, '../../src')  # noqa
import fragment_CDRs
from biopandas.pdb import PandasPdb


class TestExtractCDR(unittest.TestCase):

    # def setUp(self):

    #     # TODO put into try catch
    #     with open('test_config.yaml', 'r') as yaml_file:
    #         config_data = yaml.safe_load(yaml_file)

    def test_get_CDR_frag_dict(self):
        
        fragment_dict = fragment_CDRs.get_CDR_frag_dict()
        actual_frag_seq = 'SQGV'
        self.assertEqual(actual_frag_seq, fragment_dict['sequence'][0])

    def test_get_config_data(self):

        dirs = fragment_CDRs.get_config_data(path_to_yaml = '../../run/config.yaml')
        self.assertIsNotNone(dirs)

    
if __name__ == '__main__':
    unittest.main()
