#!/usr/bin/env python

from Bio.PDB import *
import numpy
import sys
import glob
import itertools

def replaceRes(temp_struct, temp_chain, temp_res_id, opt_struct, opt_chain, opt_res_id):
  opt_res = opt_struct[0][opt_chain][opt_res_id]
  insert_pos = [res for res in temp_struct[0][temp_chain]].index(temp_struct[0][temp_chain][temp_res_id])
  temp_struct[0][temp_chain].detach_child((' ', temp_res_id, ' '))
  temp_struct[0][temp_chain].insert(insert_pos, opt_res)
  return temp_struct
  
if __name__ == "__main__":
  parser = PDBParser(QUIET=True)
  io = PDBIO()

  qdir = sys.argv[1]
  cdr_code = qdir.split("/")[-1].split("_EpiAb")[0]
  epiAb_struct = parser.get_structure('epiAb', qdir+'/'+cdr_code+'_EpiAb.pdb')
  io.set_structure(epiAb_struct[0]["B"])
  io.save(qdir+'/'+cdr_code+'.cdr.pdb')
  
  opt_dict = {}
  for opt_qdir in glob.glob(qdir+'/*contact'):
    ab_opt = int(opt_qdir.split("contact")[0].split("_")[-1])
    opt_dict[ab_opt] = []
    for opt_path in glob.glob(opt_qdir+'/*_optimized.pdb'):
      mut_pos = opt_path.split("/")[-1].split("_")[2]
      try:
        if type(int(mut_pos)) is int:
          pass
      except:
        opt_dict[ab_opt].append(opt_path)
  opt_dict = dict(sorted({k:v for k,v in opt_dict.items() if v}.items()))
  opt_list = list()
  for idx, a in enumerate(opt_dict):
    opt_list.extend(opt_dict[a])
    for b in list(opt_dict.keys())[idx + 1:]:
      opt_list.extend(list(itertools.product(opt_dict[a], opt_dict[b])))

  ab_struct = parser.get_structure('original', qdir+'/'+cdr_code+'.cdr.pdb')
  for opt_path in opt_list:
    temp_struct = ab_struct.copy()
    if type(opt_path) is str:
      mut_pos = opt_path.split("/")[-1].split("_")[2]
      mut_pos_id = int(opt_path.split("/")[-1].split("_")[2][1:-1])
      opt_struct = parser.get_structure('original', opt_path)
      new_opt_strct = replaceRes(temp_struct, "B", mut_pos_id, opt_struct, "B", mut_pos_id)
      io.set_structure(temp_struct[0]["B"])
      io.save(qdir+'/'+cdr_code+'_'+mut_pos+'.cdr.pdb')
      print(qdir+'/'+cdr_code+'_'+mut_pos+'.cdr.pdb')
    else:
      sum_mut_pos = []
      for c, sub_opt_path in enumerate(opt_path):
        mut_pos = sub_opt_path.split("/")[-1].split("_")[2]
        sum_mut_pos.append(mut_pos)
        mut_pos_id = int(sub_opt_path.split("/")[-1].split("_")[2][1:-1])
        opt_struct = parser.get_structure('original', sub_opt_path)
        new_opt_strct = replaceRes(temp_struct, "B", mut_pos_id, opt_struct, "B", mut_pos_id)
        temp_struct = new_opt_strct.copy()
      io.set_structure(new_opt_strct[0]["B"])
      io.save(qdir+'/'+cdr_code+'_'+str('-'.join(sorted(sum_mut_pos)))+'.cdr.pdb')
      print(qdir+'/'+cdr_code+'_'+str('-'.join(sum_mut_pos))+'.cdr.pdb')
