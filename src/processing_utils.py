import os


def checkDir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    return dir_path


def dropDup(x):
  return list(dict.fromkeys(x))


def get_list_contains_str(lst, string):
    return [val for val in lst if string in val]


def call_d3to1():
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    return d3to1


def get_d3to1(aa_list):
    seq = []
    d3to1 = call_d3to1()
    for aa in aa_list:
        seq.append(d3to1[aa])
    return ''.join(seq)


def get_lineinlist(file):
    with open(file, 'r') as f:
        lines = [line.rstrip() for line in f]
    return lines
