import pickle


def open_pickle(filename):
    """
    Open a pickle file and load the data.

    Parameters:
    - filename (str): The name of the pickle file to be opened.

    Returns:
    - data: The deserialized data from the pickle file.

    Raises:
    - FileNotFoundError: If the specified file does not exist.
    - IOError: An error occurs while reading the file or unpickling the data.
    - pickle.UnpicklingError: If there is an issue with unpickling the data.
    """

    try:
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        return data
    except FileNotFoundError:
        raise FileNotFoundError(f"The file '{filename}' does not exist.")
    except IOError as e:
        raise IOError(f"Error reading the file '{filename}': {e}")
    except pickle.UnpicklingError as e:
        raise pickle.UnpicklingError(f"Error unpickling data from "
                                     f"'{filename}': {e}")


def refactor_abag(ids_list, ab_list, ag_list):
    """
    Combine lists of PDB IDs, AB-like regions, and AG-like regions into a dict

    Parameters:
    - ids_list (list): List of PDB IDs.
    - ab_list (list): List of AB-like regions.
    - ag_list (list): List of AG-like regions.

    Returns:
    - dict: A dictionary with PDB IDs as keys and corresponding AB-AG pairs as
     a list of tuples.

    Raises:
    - ValueError: If the lengths of the input lists are not equal.
    """

    if len(ids_list) != len(ab_list) or len(ids_list) != len(ag_list):
        raise ValueError("Input lists must have the same length.")

    abag_dict = {}

    # Iterate through the lists
    for id_, current_id in enumerate(ids_list):
        current_ab = ab_list[id_]
        current_ag = ag_list[id_]

        # Check if the ID is already in the dictionary
        if current_id in abag_dict:
            # Append the tuple to the existing list
            abag_dict[current_id].append((current_ab, current_ag))
        else:
            # Create a new list with the tuple for the new ID
            abag_dict[current_id] = [(current_ab, current_ag)]

    return abag_dict


def parse_string_to_list_ints(input_string):
    """
    Parse a string representation of regions into a list of position indices.

    Parameters:
    - input_string (str): Regions represented as strings with the format:
      '300|A|TRP,301|A|LEU,302|A|PRO,303|A|LEU,304|A|GLY'
      'position_index | chain | amino acid, position_index | chain ...'

    Returns:
    - index_list: List of position indices as integers.
    For example returns [300, 301, 302, 303, 304] for the string above.

    Raises:
    - ValueError: If the input string is not in the expected format.
    """

    try:
        # Split the input string by ","
        items = input_string.split(',')

        # Extract the first part of each item and convert it to an integer
        index_list = [int(item.split('|')[0]) for item in items]
        return index_list
    except (ValueError, IndexError):
        raise ValueError("Invalid input string format. Expected format: "
                         "'position_index|chain|amino acid name,...'")


def parse_string_to_list_chain_letters(input_string):
    """
    Parse a string representation of regions into a list of chain letters.

    Parameters:
    - input_string (str): Regions represented as strings with the format:
      '300|A|TRP,301|A|LEU,302|A|PRO,303|A|LEU,304|A|GLY'
      'position_index | chain | amino acid, position_index | chain ...'

    Returns:
    - chain_list: list of chain letters.
    This would be ['A','A','A','A','A'] for the example string above

    Raises:
    - ValueError: If the input string is not in the expected format.
    """

    try:
        # Split the input string by ","
        items = input_string.split(',')

        # Extract the second part of each item (chain) and create a list
        chain_list = [item.split('|')[1] for item in items]
        return chain_list
    except (ValueError, IndexError):
        raise ValueError("Invalid input string format. Expected format: "
                         "'position_index|chain|amino acid name,...'")

