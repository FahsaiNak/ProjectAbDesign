import sys
sys.path.insert(0, "../../src/")  # noqa
from PDB90 import download_and_uncompress_files, rename_files
import unittest
import os
from unittest import mock
from Bio.PDB import PDBParser, PDBIO


class TestPDB90(unittest.TestCase):

    @mock.patch('os.makedirs')
    @mock.patch('subprocess.run')
    def test_folder_creation(self, mock_run, mock_makedirs):
        # Test if the function creates the folder when it doesn't exist
        download_and_uncompress_files('test_folder', 'test.csv')
        mock_makedirs.assert_called_once_with('test_folder', exist_ok=True)

    @mock.patch('os.makedirs')
    @mock.patch('subprocess.run')
    def test_download_and_uncompress_files(self, mock_run, mock_makedirs):
        # Test if the files are being correctly downloaded and uncompressed
        with open('test.csv', 'w') as f:
            f.write('1234,test\n')
        download_and_uncompress_files('test_folder', 'test.csv')
        mock_run.assert_any_call(["wget",
                                  "-P",
                                  'test_folder',
                                  "ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb1234.ent.gz"])  # noqa
        mock_run.assert_any_call(["gunzip",
                                  os.path.join('test_folder',
                                               "pdb1234.ent.gz")])

    @mock.patch('glob.glob')
    @mock.patch('os.rename')
    def test_rename_files(self, mock_rename, mock_glob):
        # Test if the files are being correctly renamed
        mock_glob.return_value = ['test_folder/test1.ent',
                                  'test_folder/test2.ent']
        rename_files('test_folder')
        calls = [mock.call('test_folder/test1.ent', 'test_folder/t1.e.pdb'),
                 mock.call('test_folder/test2.ent', 'test_folder/t2.e.pdb')]
        mock_rename.assert_has_calls(calls)


if __name__ == '__main__':
    unittest.main()
