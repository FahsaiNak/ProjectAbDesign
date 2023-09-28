#!/usr/bin/env python

import os
import numpy
from Bio.PDB import *
import sys

class ResSelect(Select):
    def __init__(self,select_res):
        self.select_res = select_res
    def accept_residue(self, res):
        if res in self.select_res:
            return 1
        else:
            return 0

if __name__ == "__main__":
  parser = PDBParser(QUIET=True)
  io = PDBIO()
  
  full_path = sys.argv[1]
  wk_dir = sys.argv[2]
  full_struct = parser.get_structure('EpiAb', full_path)
  pep = list(full_struct[0]["A"].get_residues())[180:]
  for res in range(len(pep)-3):
    window_size = 4+res
    #print('Fragment size:',window_size)
    for i in range(len(pep) - window_size+1):
      element = list(pep)[i: i + window_size]
      #print(len(element), element)
      save_dir = wk_dir+'Ag'+'_res'+str(window_size)+'_'+str(element[0].id[1])+'-'+str(element[-1].id[1])
      os.mkdir(save_dir)
      io.set_structure(full_struct[0]["A"])
      io.save(save_dir+'/'+'Ag'+'_res'+str(window_size)+'_'+str(element[0].id[1])+'-'+str(element[-1].id[1])+'.pdb', ResSelect(element))