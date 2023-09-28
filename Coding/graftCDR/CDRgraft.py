#!/usr/bin/env python

from Bio.PDB import *
import numpy
import sys
import os

def replaceRes(temp_struct, temp_chain, temp_res, opt_struct, opt_chain, opt_res):
  int_res = opt_res.copy()
  insert_pos = [res for res in temp_struct[0][temp_chain]].index(temp_res)
  temp_struct[0][temp_chain].detach_child(temp_res.id)
  int_res.id = tuple(temp_res.id)
  temp_struct[0][temp_chain].insert(insert_pos, int_res)
  
if __name__ == "__main__":
  parser = PDBParser(QUIET=True)
  io = PDBIO()
  
  match_path = sys.argv[1]
  
  cdr_path = sys.argv[2]
  cdr_name = cdr_path.split("/")[-1].split(".")[0]
  cdr = parser.get_structure('cdr', cdr_path)
  cdr_chain = [chain.id for chain in cdr[0].get_chains()][0]
  
  TCRtype = sys.argv[3]
  save_dir = sys.argv[4]
  
  matchfile = match_path.split("/")[-1]
  match_struc = parser.get_structure('match', match_path)
  match_chain = [chain.id for chain in match_struc[0].get_chains()][0]
  full = parser.get_structure('full', "./TCRmatch/"+cdr_name+"/"+TCRtype+"_lowcutoff.match.full.struct/"+matchfile)
  full_chain = [chain.id for chain in full[0].get_chains()][0]
  #print([res for res in full[0][[chain.id for chain in full[0].get_chains()][0]]])
  for match_res, cdr_res in zip(match_struc[0][match_chain], cdr[0][cdr_chain]):
    replaceRes(full, full_chain, match_res, cdr, cdr_chain, cdr_res)
  #print([res for res in full[0][[chain.id for chain in full[0].get_chains()][0]]])
  io.set_structure(full[0])
  io.save(save_dir+"/"+cdr_name+".pdb")

  match_nres_list = list()
  for res in match_struc[0][match_chain]:
    if res.id[2] != " ":
      match_nres_list.append(str(res.id[1])+res.id[2])
    else:
      match_nres_list.append(str(res.id[1]))
  report_file = save_dir+"/graft_report.txt"
  if os.path.exists(report_file):
    with open(report_file, "a") as File:
      File.write(cdr_name+"\t"+str(match_nres_list)+'\n')
  else:
    with open(report_file, "w") as File:
      File.write(cdr_name+"\t"+str(match_nres_list)+'\n')
    
  
  