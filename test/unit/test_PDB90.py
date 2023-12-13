import os
import shutil
import unittest
import sys
import gzip
sys.path.insert(0, "../../src/")  # noqa
from PDB90 import download_files, uncompress_files, rename_files


class TestScript(unittest.TestCase):
    """Test cases for the script."""

    def setUp(self):
        """Create a test directory."""
        os.makedirs('test_folder', exist_ok=True)

    def tearDown(self):
        """Clean up the test directory."""
        shutil.rmtree('test_folder')

    def test_folder_creation(self):
        """Test if the function creates the folder when it doesn't exist."""
        download_files('test_folder', '../Datasets/test.csv')
        self.assertTrue(os.path.exists('test_folder'))

    def test_download_files(self):
        """Test if the files are being correctly downloaded."""
        with open('../Datasets/test.csv', 'w') as f:
            f.write('6anr,test\n')  # Use a valid PDB ID for testing
        download_files('test_folder', '../Datasets/test.csv')
        # Check if the .gz file exists
        self.assertTrue(any(os.path.isfile(os.path.join('test_folder', f))
                            and f.endswith('.gz') for f in os.listdir('test_folder')))  # noqa

    def test_uncompress_files(self):
        """Test if the files are being correctly uncompressed."""
        # Create some test files
        content = b'Test content'
        with gzip.open('test_folder/pdb6anr.ent.gz', 'wb') as f:
            f.write(content)
        uncompress_files('test_folder')
        # Check if the files are uncompressed
        self.assertTrue(os.path.exists('test_folder/pdb6anr.ent'))

    def test_rename_files(self):
        """Test if the files are being correctly renamed."""
        # Create some test files
        open('test_folder/pdb6anr.ent', 'a').close()
        rename_files('test_folder')
        # Check if the files are renamed
        self.assertTrue(os.path.exists('test_folder/6anr.pdb'))


if __name__ == '__main__':
    unittest.main()
