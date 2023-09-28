#!/usr/bin/env python

import numpy
from Bio.PDB import *
from Bio.PDB.SASA import ShrakeRupley
import sys
import os

def replaceRes(temp_struct, temp_chain, temp_res, opt_struct, opt_chain, opt_res):
  int_res = opt_res.copy()
  insert_pos = [res for res in temp_struct[0][temp_chain]].index(temp_res)
  temp_struct[0][temp_chain].detach_child(temp_res.id)
  int_res.id = tuple(temp_res.id)
  temp_struct[0][temp_chain].insert(insert_pos, int_res)

def lines_that_contain(string, fp):
    return [line for line in fp if string in line]

def findIndex(mark_list):
  index_list = list()
  for i,mark in enumerate(mark_list):
    try:
      if type(int(mark)) == int:
        index_list.append(i)
    except:
      pass
  return index_list[0], index_list[-1]

def extractResID(graft_report_file, name):
  cdr_res_list = []
  with open(graft_report_file, "r") as fp:
    for line in lines_that_contain(name, fp):
      for ele in line.strip().split("\t")[1].split(", "):
        firstind, lastind = findIndex([_ for _ in ele])
        if ele[lastind+1] != "'":
          cdr_res_list.append(tuple((' ', int(ele[firstind:lastind+1]), ele[lastind+1])))
        else:
          cdr_res_list.append(tuple((' ', int(ele[firstind:lastind+1]), ' ')))
  return cdr_res_list

if __name__ == "__main__":
  parser = PDBParser(QUIET=True)
  io = PDBIO()
  sup = Superimposer()
  
  cdr_name = sys.argv[1]
  neo_resi = int(sys.argv[2])
  wk_dir = sys.argv[3]
  template_name = wk_dir.split("/")[-1]
  template_pdb = wk_dir.split("/")[-1].split("_")[0]
  template_b_chain = wk_dir.split("/")[-1].split("_")[1]
  #print(template_name)
  
  #compute SASA from the target residues before coupling with CDR-like-TCRvb
  target_struct = parser.get_structure('target', wk_dir+"/target_aligned.pdb")
  target_chain = list(target_struct[0].get_chains())[0].id
  sr = ShrakeRupley()
  sr.compute(target_struct[0], level="R")
  SASA_store = dict()
  for res in target_struct[0][target_chain]:
    SASA_store[res.id[1]] = round(res.sasa,2)
  #print(SASA_store)
  
  #compute SASA from the target residues after coupling with template TCRvb
  pmhctcrvb_struct = parser.get_structure('TCRvb', wk_dir+"/target_pMHCTCRvb.pdb")
  sr.compute(pmhctcrvb_struct[0], level="R")
  decreasedSASA_nres_control = []
  for res in pmhctcrvb_struct[0][target_chain]:
    if SASA_store[res.id[1]] >= round(res.sasa,2) + 0.4:
      decreasedSASA_nres_control.append(res.id[1])
  #print(decreasedSASA_nres_control)

  #graft the CDR-like-TCRb3 region to the single variable domain from a TCR beta chain template
  template_vb_struct = parser.get_structure('TCRvb', wk_dir+"/"+template_name+"_TCRvb.pdb")
  temp_struct = template_vb_struct.copy()
  graft_struct = parser.get_structure('TCRvb', wk_dir+"/"+cdr_name+".pdb")
  graft_nres_list = [res.id for res in graft_struct[0][template_b_chain]]
  ref_atoms = []
  for ref_res in temp_struct[0][template_b_chain]:
    if ref_res.id in graft_nres_list:
      ref_atoms.append(ref_res['CA'])
  alt_atoms = []
  for alt_res in graft_struct[0][template_b_chain]:
    alt_atoms.append(alt_res['CA'])
  sup.set_atoms(ref_atoms, alt_atoms)
  sup.apply(graft_struct[0].get_atoms())   
  for res in temp_struct[0][template_b_chain]:
    if res.id in graft_nres_list:
      replaceRes(temp_struct, template_b_chain, res, graft_struct, template_b_chain, graft_struct[0][template_b_chain][res.id])
  io.set_structure(temp_struct[0][template_b_chain])
  io.save(wk_dir+"/TCRb3-"+cdr_name+".pdb")
  
  #combine with target pMHC
  target_struct_temp = target_struct.copy()
  grafted_tcrvb_struct = parser.get_structure('graftedTCRvb', wk_dir+"/TCRb3-"+cdr_name+".pdb")
  target_struct_temp[0].add(grafted_tcrvb_struct[0][template_b_chain])
  io.set_structure(target_struct_temp[0])
  io.save(wk_dir+"/TCRb3-"+cdr_name+"_pMHCTCRvb.pdb")
  
  #compute SASA from the target residues after coupling with CDR-like-TCRvb
  complex_struct = parser.get_structure('TCRvb', wk_dir+"/TCRb3-"+cdr_name+"_pMHCTCRvb.pdb")
  sr.compute(complex_struct[0], level="R")
  SASA_complex_store = dict()
  decreasedSASA_nres = []
  for res in complex_struct[0][target_chain]:
    SASA_complex_store[res.id[1]] = round(res.sasa,2)
    if SASA_store[res.id[1]] >= round(res.sasa,2) + 0.4:
      decreasedSASA_nres.append(res.id[1])
  #print(decreasedSASA_nres)
  
  #compute SASA from the target residues after detaching cdr-like residues
  neoint = False
  cdr_nres_list = extractResID(wk_dir+"/graft_report.txt", cdr_name)
  print(cdr_name, cdr_nres_list)
  crop_struct = parser.get_structure('TCRvb', wk_dir+"/TCRb3-"+cdr_name+"_pMHCTCRvb.pdb")
  for res in crop_struct[0][template_b_chain]:
    if res.id in cdr_nres_list:
      crop_struct[0][template_b_chain].detach_child(res.id)
  increasedSASA_nres = []
  sr.compute(crop_struct[0], level="R")
  for res in crop_struct[0][target_chain]:
    if round(res.sasa,2) >= SASA_complex_store[res.id[1]]+0.4:
      increasedSASA_nres.append(res.id[1])
  if len(increasedSASA_nres) > 0:
    print("pMHC-CDR contact:",increasedSASA_nres)
    if neo_resi in increasedSASA_nres:
      neoint = True
      
  print(template_name+"\tTCRb3-"+cdr_name+"\t"+str(len(decreasedSASA_nres))+"\t"+str(len(decreasedSASA_nres_control))+"\t"+str(len(set(decreasedSASA_nres_control)&set(decreasedSASA_nres)))+"\t"+str(list(set([n for n in range(181,190)])&set(decreasedSASA_nres))), neoint)
  if os.path.exists(wk_dir+"/interaction_report.txt"):
    with open(wk_dir+"/interaction_report.txt", "a") as File:
      File.write("TCRb3-"+cdr_name+"\t"+template_name+"\t"+str(len(set(decreasedSASA_nres_control)&set(decreasedSASA_nres)))+"\t"+str(len(set(decreasedSASA_nres_control)-set(decreasedSASA_nres)))+"\t"+str(neoint)+'\n')
  else:
    with open(wk_dir+"/interaction_report.txt", "w") as File:
      File.write("TCRb3-"+cdr_name+"\t"+template_name+"\t"+str(len(set(decreasedSASA_nres_control)&set(decreasedSASA_nres)))+"\t"+str(len(set(decreasedSASA_nres_control)-set(decreasedSASA_nres)))+"\t"+str(neoint)+'\n')
