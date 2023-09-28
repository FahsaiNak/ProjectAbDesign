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
  TCR = sys.argv[2]
  tcrm_name, pdb_template  = line.split()[0], line.split()[1]
  cdr_name = tcrm_name[6:]
  data_dir = "./TCRgraft/"+TCR+"/"+pdb_template
  
  #compute SASA from the peptide residues when coupling with CDR-grafted-TCRm
  complex_struc = parser.get_structure('pMHCTCRbv', data_dir+"/"+tcrm_name+"_pMHCTCRvb.pdb")
  pep_chain = list(complex_struc[0].get_chains())[0].id
  tcr_chain = list(complex_struc[0].get_chains())[1].id
  sr.compute(complex_struc[0], level="R")
  SASA_store_pep = dict()
  for res in complex_struc[0][pep_chain]:
    SASA_store_pep[res.id[1]] = round(res.sasa,2)
  SASA_store_tcr = dict()
  for res in complex_struc[0][tcr_chain]:
    SASA_store_tcr[res.id[1]] = round(res.sasa,2)
  
  tcr_struc = parser.get_structure('MHCTCRvb', data_dir+"/"+tcrm_name+"_pMHCTCRvb.pdb")
  for res in complex_struc[0][pep_chain]:
    if res.id[1] in range(181, 190): #9-length peptide
      tcr_struc[0][pep_chain].detach_child(res.id)
  sr.compute(tcr_struc[0], level="R")
  res_contact_cdr = list()
  for res in tcr_struc[0][tcr_chain]:
    if res.id[1] in [_ for _ in range(105,118)]: #CDR3 IMGT-numbering
      if round(res.sasa,2) >= SASA_store_tcr[res.id[1]]+0.4:
        res_contact_cdr.append(res)
  for cdr_res in res_contact_cdr:
    print(cdr_res)
    cdr_res_opt = str(cdr_res.id[1])
    if cdr_res.id[2] != " ":
      cdr_res_opt = str(cdr_res.id[1])+cdr_res.id[2]
    crop_struc = parser.get_structure('crop_pMHCTCRbv', data_dir+"/"+tcrm_name+"_pMHCTCRvb.pdb")
    crop_struc[0][tcr_chain].detach_child(cdr_res.id)
    sr.compute(crop_struc[0], level="R")
    res_contact_target = list()
    for pep_res in crop_struc[0][pep_chain]:
      if pep_res.id[1] in range(181, 190):
        if round(pep_res.sasa,2) >= SASA_store_pep[pep_res.id[1]]+0.4:
          res_contact_target.append(pep_res)
    if len(res_contact_target) > 0:
      selected_res = res_contact_target+res_contact_cdr
      print(selected_res)
      io.set_structure(complex_struc[0])
      os.mkdir("./TCRopt/"+TCR+"/"+pdb_template+"/"+tcrm_name+"/"+tcrm_name+"_"+cdr_res_opt+"contact")
      io.save("./TCRopt/"+TCR+"/"+pdb_template+"/"+tcrm_name+"/"+tcrm_name+"_"+cdr_res_opt+"contact/"+tcrm_name+"_"+cdr_res_opt+"contact.pdb", ResSelect(selected_res))
    print("--------")