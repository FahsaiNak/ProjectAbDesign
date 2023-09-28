#!/usr/bin/env python

import os
import numpy
from Bio.PDB import *
from Bio.PDB.SASA import ShrakeRupley
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
  sr = ShrakeRupley()
  
  line = sys.argv[1]
  ind, pdb, Ab_chain = str(line.split()[0]), str(line.split()[1]), str(line.split()[2])
  Ab_pos = [int(n) for n in line.split("[")[1].split("] ")[0].split(", ")]
  Ag_chain = line.split("[")[1].split("] ")[1].split()[0]
  Ag_pos = [int(n) for n in line.split("[")[2].split("] ")[0].split(", ")]
  epi = line.split("[")[2].split("] ")[1].split()[0]
  epi_pos = [int(n) for n in line.split("[")[3].split("]")[0].split(", ")]
  #print(ind, pdb, Ab_chain, Ab_pos, Ag_chain, Ag_pos, epi, epi_pos)
  Name = ind+"_Ag"+str(Ag_pos[0])+"-"+str(Ag_pos[-1])+"_EpiAb"
  epi_chain = "A"
  Ab_chain = "B"
  
  #compute SASA from the epitope residues when coupling with full CDR-like
  full_struc = parser.get_structure('EpiAb', epi+"/"+Name+"/"+Name+".pdb")
  sr.compute(full_struc[0], level="R")
  SASA_store = dict()
  for res in full_struc[0][epi_chain]:
    SASA_store[res.id[1]] = round(res.sasa,2)
  
  #compute SASA from the target residues after detaching each cdr-like residue
  for Ab_res in full_struc[0][Ab_chain]:
    crop_struc = parser.get_structure('EpiAb', epi+"/"+Name+"/"+Name+".pdb")
    crop_struc[0][Ab_chain].detach_child(Ab_res.id)
    increasedSASA_nres = []
    sr.compute(crop_struc[0], level="R")
    for res in crop_struc[0][epi_chain]:
      if round(res.sasa,2) >= SASA_store[res.id[1]]+0.4:
        increasedSASA_nres.append(res)
    if len(increasedSASA_nres) > 0:
      print([res.id for res in increasedSASA_nres])
      selected_res = increasedSASA_nres+[res for res in full_struc[0][Ab_chain]]
      io.set_structure(full_struc[0])
      os.mkdir(epi+"/"+Name+"/"+Name+"_"+str(Ab_res.id[1])+"contact")
      io.save(epi+"/"+Name+"/"+Name+"_"+str(Ab_res.id[1])+"contact/"+Name+"_"+str(Ab_res.id[1])+"contact.pdb", ResSelect(selected_res))
