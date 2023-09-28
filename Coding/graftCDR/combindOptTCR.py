#!/usr/bin/env python

from Bio.PDB import *
from Bio.PDB.SASA import ShrakeRupley
import numpy
import sys
import glob
import os
import itertools
from multiprocessing import Process

def replaceRes(temp_struct, temp_chain, temp_res_id, opt_struct, opt_chain, opt_res_id):
  opt_res = opt_struct[0][opt_chain][opt_res_id]
  insert_pos = [res for res in temp_struct[0][temp_chain]].index(temp_struct[0][temp_chain][temp_res_id])
  temp_struct[0][temp_chain].detach_child(temp_res_id)
  temp_struct[0][temp_chain].insert(insert_pos, opt_res)
  return temp_struct

def combinedOptTCR(tcr_struct, tcr_chain, tcr_code, com, qdir, pmhc_struct, SASA_pmhc, BSA_ogcom, pdb):
  temp_struct = tcr_struct.copy()
  sum_mut_pos = []
  for mut_pos, opt_path in com.items():
    if "contact" not in opt_path.split("/")[-1].split("_")[2]:
      sum_mut_pos.append(opt_path.split("/")[-1].split("_")[2])
    try:
      mut_pos_id = tuple((' ', int(mut_pos), ' '))
    except:
      mut_pos_id = tuple((' ', int(mut_pos[:-1]), mut_pos[-1]))
    #print(mut_pos,mut_pos_id,opt_path)
    opt_struct = parser.get_structure('opt', opt_path)
    new_opt_strct = replaceRes(temp_struct, tcr_chain, mut_pos_id, opt_struct, tcr_chain, mut_pos_id)
    temp_struct = new_opt_strct.copy()
  io.set_structure(temp_struct[0][tcr_chain])
  if len(sum_mut_pos) > 0:
    opt_tcr_code = tcr_code+'_'+str('-'.join(sorted(sum_mut_pos)))
  else:
    opt_tcr_code = tcr_code
  io.save(qdir+'/'+opt_tcr_code+'.pdb')
  
  opt_tcr = parser.get_structure('optimized', qdir+'/'+opt_tcr_code+'.pdb')
  sr.compute(opt_tcr, level="S")
  SASA_tcr = opt_tcr.sasa
  complex_struct = pmhc_struct.copy()
  complex_struct[0].add(opt_tcr[0][tcr_chain])
  sr.compute(complex_struct, level="S")
  SASA_com = complex_struct.sasa
  BSA_com = (SASA_pmhc+SASA_tcr)-SASA_com
  #print(opt_tcr_code+"-SASA", round(SASA_tcr,2), "\nComplex-SASA", round(SASA_com,2), "\nBSA", round(BSA_com,2))
  if round(BSA_com,2) < round(BSA_ogcom,2):
    os.remove(qdir+'/'+opt_tcr_code+'.pdb')

if __name__ == "__main__":
  parser = PDBParser(QUIET=True)
  io = PDBIO()
  sr = ShrakeRupley()

  qdir = sys.argv[1]
  tcr_code = qdir.split("/")[-1]
  tcrdir = "./TCRgraft/"+"/".join(qdir.split("/")[2:-1])
  pdb = qdir.split("/")[3]
  tcr_chain = pdb.split("_")[-1]
  
  opt_dict = {}
  for opt_qdir in glob.glob(qdir+'/*contact'):
    ab_opt = opt_qdir.split("contact")[0].split("_")[-1]
    opt_dict[ab_opt] = [opt_qdir+"/"+opt_qdir.split("/")[-1]+".pdb"]
    for opt_path in glob.glob(opt_qdir+'/*_optimized.pdb'):
      mut_pos = opt_path.split("/")[-1].split("_")[2]
      try:
        if type(int(mut_pos[0])) is int:
          pass
      except:
        opt_dict[ab_opt].append(opt_path)
  mut_dict = [{k: v for k, v in zip(sorted(opt_dict.keys()), list_prod_value)} for list_prod_value in itertools.product(*(opt_dict[k] for k in sorted(opt_dict.keys())))]
  
  tcr_struct = parser.get_structure('original', tcrdir+'/'+tcr_code+'.pdb')
  sr.compute(tcr_struct, level="S")
  SASA_ogtcr = tcr_struct.sasa
  
  pmhc_struct = parser.get_structure('pmhc', tcrdir+'/target_aligned.pdb')
  sr.compute(pmhc_struct, level="S")
  SASA_pmhc = pmhc_struct.sasa 
  
  complex_struct = pmhc_struct.copy()
  complex_struct[0].add(tcr_struct[0][tcr_chain])
  sr.compute(complex_struct, level="S")
  SASA_ogcom = complex_struct.sasa
  BSA_ogcom = (SASA_pmhc+SASA_ogtcr)-SASA_ogcom
  #print("pMHC-SASA", round(SASA_pmhc,2), "\nogTCR-SASA", round(SASA_ogtcr,2), "\nComplex-SASA", round(SASA_ogcom,2), "\nBSA", round(BSA_ogcom,2))
  for com in mut_dict:
    process = Process(target=combinedOptTCR, args=(tcr_struct, tcr_chain, tcr_code, com, qdir, pmhc_struct, SASA_pmhc, BSA_ogcom, pdb))
    process.start()
