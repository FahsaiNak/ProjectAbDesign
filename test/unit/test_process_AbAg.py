import unittest
import sys
import pickle
import pandas as pd
sys.path.insert(0, '../../src')  # noqa
import process_AbAg_utils as ut


class TestRefactor(unittest.TestCase):

    def test_open_pickle_success(self):
        # test_pickle in data folder, contains test_data
        test_filename = \
            '../Datasets/test_pkl_and_regions_for_vis/test_pickle.pkl'
        test_data = {'example_key': 'example_value'}

        try:
            # Test if the function successfully opens and loads the pickle file
            loaded_data = ut.open_pickle(test_filename)
            self.assertEqual(loaded_data, test_data)
        except Exception as e:
            self.fail(f"Unexpected exception: {e}")

    def test_open_pickle_file_not_found(self):
        # Test if the function raises FileNotFoundError for a non-existing file
        with self.assertRaises(FileNotFoundError):
            ut.open_pickle('non_existing_file.pkl')

    def test_open_pickle_unpickling_error(self):
        # Create an empty file (not a valid pickle file)
        test_filename = \
            '../Datasets/test_pkl_and_regions_for_vis/corrupted_pickle.pkl'

        # Test if the function raises UnpicklingError for an invalid pickle
        # file
        with self.assertRaises(pickle.UnpicklingError):
            ut.open_pickle(test_filename)

    def test_parse_string_to_list_ints_success(self):
        # Test if the function successfully parses the input string
        input_string = '300|A|TRP,301|A|LEU,302|A|PRO,303|A|LEU,304|A|GLY'
        expected_result = [300, 301, 302, 303, 304]
        result = ut.parse_string_to_list_ints(input_string)
        self.assertEqual(result, expected_result)

    def test_parse_string_to_list_ints_invalid_format(self):
        # Test if the function raises ValueError for an invalid input string
        # format
        with self.assertRaises(ValueError):
            ut.parse_string_to_list_ints('invalid_format_string')

    def test_parse_string_to_list_chain_letters_success(self):
        # Test if the function successfully parses the input string
        input_string = '300|A|TRP,301|A|LEU,302|A|PRO,303|A|LEU,304|A|GLY'
        expected_result = ['A', 'A', 'A', 'A', 'A']
        result = ut.parse_string_to_list_chain_letters(input_string)
        self.assertEqual(result, expected_result)

    def test_parse_string_to_list_invalid_format(self):
        # Test if the function raises ValueError for an invalid input string
        # format
        with self.assertRaises(ValueError):
            ut.parse_string_to_list_chain_letters('invalid_format_string')

    def test_filter_df_success(self):
        """
        Test filtering DataFrame with existing PDB-IDs in the CSV file.
        """
        unfiltered_df = pd.DataFrame({
            'pdb_id': ['1abc', '1abc', '2xyz', '3def'],
            'value': [10, 47, 20, 30]
        })
        test_filename = \
            '../Datasets/test_pkl_and_regions_for_vis/test_filter_pdb_id.csv'
        # function call
        filtered_df = ut.filter_df(unfiltered_df, test_filename)
        self.assertEqual(len(filtered_df), 2)

    def test_filter_df_file_not_found(self):
        """
        Test FileNotFoundError when CSV file is not found.
        """
        unfiltered_df = pd.DataFrame({
            'pdb_id': ['1abc', '1abc', '2xyz', '3def'],
            'value': [10, 47, 20, 30]
        })
        wrong_filename = 'not_a_thing.csv'

        # Call the filter_df function with a non-existing file
        with self.assertRaises(FileNotFoundError):
            ut.filter_df(unfiltered_df, wrong_filename)

    def test_filter_df_no_matches_found(self):
        """
        Test ValueError when no matching PDB-IDs are found.
        """
        unfiltered_df = pd.DataFrame({
            'pdb_id': ['2abc', '2abc', '2xyz', '3def'],
            'value': [10, 47, 20, 30]
        })
        test_filename = \
            '../Datasets/test_pkl_and_regions_for_vis/test_filter_pdb_id.csv'

        with self.assertRaises(ValueError):
            ut.filter_df(unfiltered_df, test_filename)


if __name__ == '__main__':
    unittest.main()
