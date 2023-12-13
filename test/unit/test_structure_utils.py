import sys
import unittest
import random
import os
from Bio.PDB import PDBList, PDBParser, ShrakeRupley
sys.path.insert(0, '../../src')
import structure_utils as stut  # noqa


class TestStrucUt(unittest.TestCase):
    def setUp(self):
        pdb_list = PDBList()
        parser = PDBParser(QUIET=True)
        sr = ShrakeRupley()
        pdb_file = os.listdir("../Datasets/PDB90_PDS")
        pdb_id_list = [_.split("/")[-1].split(".")[0] for _ in pdb_file]
        self.test_pdb_id = pdb_id_list[random.randint(0, len(pdb_id_list)-1)]
        pdb_list.download_pdb_files([self.test_pdb_id],
                                    pdir='.', file_format='pdb')
        self.pdb_file_name = f"pdb{self.test_pdb_id}.ent"  # noqa
        self.control_struct = parser.get_structure(self.test_pdb_id,
                                                   self.pdb_file_name)
        self.rand_res = random.choice(list(self.control_struct.get_residues()))
        self.rand_chain = self.rand_res.get_parent().id
        self.crop_struct = self.control_struct.copy()
        self.crop_struct[0][self.rand_chain].detach_child(self.rand_res.id)  # noqa
        self.d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',  # noqa
                      'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',  # noqa
                      'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',  # noqa
                      'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}  # noqa
        sr.compute(self.control_struct[0], level="R")
        sr.compute(self.crop_struct[0], level="R")
        self.no_sasa_struct = parser.get_structure(self.test_pdb_id,
                                                   self.pdb_file_name)
        self.agsasa_list = list()
        for chain in self.control_struct[0]:
            for res in chain:
                if res != self.rand_res:
                    sasa_val_full = round(res.sasa, 2)
                    sasa_val_crop = round(self.crop_struct[0][chain.id][res.id].sasa, 2)  # noqa
                    if sasa_val_crop >= sasa_val_full+0.4:
                        if res.resname in self.d3to1.keys():
                            self.agsasa_list.append(f"{res.id[1]}|{chain.id}|{res.resname}")  # noqa

    def test_get_structurefrompdb(self):
        structure = stut.get_structurefrompdb(self.test_pdb_id)
        self.assertEqual(list(structure.get_residues()),
                         list(self.control_struct.get_residues()))
        self.assertEqual(None, stut.get_structurefrompdb("errr"))

    def test_get_structurefromfile(self):
        structure = stut.get_structurefromfile(self.pdb_file_name)
        self.assertEqual(list(structure.get_residues()),
                         list(self.control_struct.get_residues()))
        self.assertEqual(None, stut.get_structurefromfile("errr"))
        os.remove(self.pdb_file_name)

    def test_find_SASAcontact(self):
        agsasa_list = stut.find_SASAcontact(self.control_struct,
                                            self.crop_struct)
        agsasa_empty_list = stut.find_SASAcontact(self.no_sasa_struct,
                                                  self.crop_struct)
        self.assertEqual(self.agsasa_list, agsasa_list)
        self.assertEqual([], agsasa_empty_list)
        os.remove(self.pdb_file_name)


if __name__ == '__main__':
    unittest.main()
