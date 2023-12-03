import unittest
import sys
import pickle
sys.path.insert(0, '../../src')  # noqa
import process_AbAg_utils as ut


class TestRefactor(unittest.TestCase):

    def test_open_pickle_success(self):
        # test_pickle in data folder, contains test_data
        test_filename = '../../data/test_pickle.pkl'
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
        test_filename = '../../data/corrupted_pickle.pkl'

        # Test if the function raises UnpicklingError for an invalid pickle
        # file
        with self.assertRaises(pickle.UnpicklingError):
            ut.open_pickle(test_filename)

    def test_parse_string_to_list_success(self):
        # Test if the function successfully parses the input string
        input_string = '300|A|TRP,301|A|LEU,302|A|PRO,303|A|LEU,304|A|GLY'
        expected_result = [300, 301, 302, 303, 304]

        result = ut.parse_string_to_list(input_string)
        self.assertEqual(result, expected_result)

    def test_parse_string_to_list_invalid_format(self):
        # Test if the function raises ValueError for an invalid input string
        # format
        with self.assertRaises(ValueError):
            ut.parse_string_to_list('invalid_format_string')


if __name__ == '__main__':
    unittest.main()



