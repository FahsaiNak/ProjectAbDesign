import sys
import os
import shutil
import subprocess
sys.path.insert(0, "../../src/")  # noqa
from PDB90 import download_and_uncompress_files, rename_files
import unittest
from Bio.PDB import PDBParser, PDBIO


class TestPDB90(unittest.TestCase):

    def setUp(self):
        # Create a test directory
        os.makedirs('test_folder', exist_ok=True)

    def tearDown(self):
        # Clean up the test directory
        shutil.rmtree('test_folder')

    def test_folder_creation(self):
        # Test if the function creates the folder when it doesn't exist
        download_and_uncompress_files('test_folder', '../Datasets/PDB90_test.csv')
        self.assertTrue(os.path.exists('test_folder'))

    def test_download_and_uncompress_files(self):
        # Test if the files are being correctly downloaded and uncompressed
        with open('../Datasets/PDB90_test.csv', 'w') as f:
            f.write('6anr,test\n')  # Use a valid PDB ID for testing
        download_and_uncompress_files('test_folder', '../Datasets/PDB90_test.csv')
        # Check if the .ent file exists
        self.assertTrue(any(os.path.isfile(os.path.join('test_folder', f)) and f.endswith('.ent') for f in os.listdir('test_folder')))

    def test_rename_files(self):
        # Test if the files are being correctly renamed
        # Create some test files
        open('test_folder/test1.ent', 'a').close()
        open('test_folder/test2.ent', 'a').close()
        rename_files('test_folder')
        # Check if the files are renamed
        self.assertTrue(os.path.exists('test_folder/t1.e.pdb'))
        self.assertTrue(os.path.exists('test_folder/t2.e.pdb'))


if __name__ == '__main__':
    unittest.main()
