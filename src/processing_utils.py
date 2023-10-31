import os


def checkDir(dir_path):
    """
    Check if a directory exists; if not, create it.

    Args:
        dir_path (str): The path of the directory to check/create.

    Returns:
        str: The input dir_path, whether it existed or was created.
    """
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    return dir_path


def dropDup(x):
    """
    Remove duplicates from a list while preserving the order.

    Args:
        x (list): The input list with potential duplicates.

    Returns:
        list: The list with duplicates removed.
    """
    return list(dict.fromkeys(x))


def get_list_contains_str(lst, string):
    """
    Filter a list to get elements that contain a specific substring.

    Args:
        lst (list): The input list of strings.
        string (str): The substring to search for in the elements.

    Returns:
        list: The list of elements that contain the specified substring.
    """
    return [val for val in lst if string in val]


def call_d3to1():
    """
    Create a dictionary for converting
    three-letter amino acid codes to one-letter codes.

    Returns:
        dict: A dictionary mapping three-letter codes to one-letter codes.
    """
    d3to1 = {
        'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
        'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
    }
    return d3to1


def get_d3to1(aa_list):
    """
    Convert a list of three-letter amino acid codes to one-letter codes.

    Args:
        aa_list (list): A list of three-letter amino acid codes.

    Returns:
        str: A string containing one-letter amino acid codes.
    """
    seq = []
    d3to1 = call_d3to1()
    for aa in aa_list:
        seq.append(d3to1[aa])
    return ''.join(seq)


def get_lineinlist(file):
    """
    Read lines from a file into a list.

    Args:
        file (str): The path to the file to read.

    Returns:
        list: A list of lines from the file.
    """
    with open(file, 'r') as f:
        lines = [line.rstrip() for line in f]
    return lines
