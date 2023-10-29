import sys
import unittest
import random
import os
from Bio.PDB import PDBParser, PDBList
sys.path.insert(0, '../../src')
import master_utils as mut  # noqa


class TestMasterUt(unittest.TestCase):
    def setUp(self):
        pdb_list = PDBList()
        parser = PDBParser(QUIET=True)
        pdb_file = os.listdir("../Datasets/PDB90_PDS")
        pdb_id_list = [_.split("/")[-1].split(".")[0] for _ in pdb_file]
        self.test_pdb_id = pdb_id_list[random.randint(0, len(pdb_id_list)-1)]
        pdb_list.download_pdb_files([self.test_pdb_id], pdir='.', file_format='pdb')
        self.pdb_file_name = f"pdb{self.test_pdb_id}.ent"
        self.control_struct = parser.get_structure(self.test_pdb_id, self.pdb_file_name)
        self.residue_list = list()
        for chain in self.control_struct[0]:
            for residue in chain:
                res = f"{residue.id[1]}|{chain.id}|{residue.resname}"
                self.residue_list.append(res)
    

    def test_get_all_residues(self):
        residue_list = mut.get_all_residues(self.pdb_file_name)
        self.assertEqual(self.residue_list, residue_list)
        self.assertEqual(None, mut.get_all_residues("errr.ent"))
        os.remove(self.pdb_file_name)


if __name__ == '__main__':
    unittest.main()
