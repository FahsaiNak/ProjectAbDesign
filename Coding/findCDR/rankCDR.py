#!/usr/bin/env python

import os
import numpy
from Bio.PDB import *
from Bio.PDB.SASA import ShrakeRupley
import sys
import glob
from multiprocessing import Process

def extractResid(listIDinText):
  resid = []
  for i,n in enumerate(listIDinText.split(", ")):
    if i == 0:
      resid.append(int(n.split("[")[1]))
    elif i == len(listIDinText.split(", "))-1:
      resid.append(int(n.split("]")[0]))
    else:
      resid.append(int(n))
  return resid

def findSASAcontact(full_path, cdr_chain, cdr_resid, ag_chain, ag_resid):
  sr = ShrakeRupley()
  full_struct = parser.get_structure('full', full_path)
  sr.compute(full_struct[0], level="R")
  SASA_store = dict()
  for ag_res in full_struct[0][ag_chain]:
    if ag_res.id[1] in ag_resid:
      SASA_store[ag_res.id[1]] = round(ag_res.sasa,2) 
  ag_contact = []
  for cdr_res in full_struct[0][cdr_chain]:
    if cdr_res.id[1] in cdr_resid:
      crop_struc = parser.get_structure('crop', full_path)
      crop_struc[0][cdr_chain].detach_child((' ', cdr_res.id[1], ' '))
      sr.compute(crop_struc[0], level="R")
      for res in crop_struc[0][ag_chain]:
        if res.id[1] in ag_resid:
          if round(res.sasa,2) >= SASA_store[res.id[1]]+0.4:
            ag_contact.append(res.id[1])
  return sorted(list(dict.fromkeys(ag_contact)))

def parallelRun(code_id, opt_path, cdr_resid, fullepitope_path, qdir, epi, converted_ag_contact_native, native_target, neo_resid):
  opt_struct = parser.get_structure('opt', opt_path)
  opt_cdr_resid = [res.id[1] for res in Selection.unfold_entities(opt_struct[0]["B"],"R") if res.id[1] in cdr_resid]
  target_struct = parser.get_structure('target', fullepitope_path)
  target_struct[0].add(opt_struct[0]["B"])
  io.set_structure(target_struct[0])
  io.save(qdir+'/pdb/'+code_id+'.EpiAb.pdb')
  ag_contact_target = findSASAcontact(qdir+'/pdb/'+code_id+'.EpiAb.pdb', "B", opt_cdr_resid, "A", [res.id[1] for res in target_struct[0]["A"].get_residues()])
  os.remove(qdir+'/pdb/'+code_id+'.EpiAb.pdb')
  ag_contact_epi = set(native_target.values())&set(ag_contact_target)
  print("Target",ag_contact_target)
  shared_int = set(converted_ag_contact_native)&set(ag_contact_target)
  print("Shared",shared_int, len(shared_int))
  nonshared_int = ag_contact_epi-set(converted_ag_contact_native)
  print("Nonshared",nonshared_int, len(nonshared_int))
  if neo_resid in ag_contact_target and neo_resid in converted_ag_contact_native:
    neo_int = True
    print("Found neo-site interaction!")
  else:
    neo_int = False
  out = opt_path.split("/")[-1]+'\t'+epi+'\t'+str(len(shared_int))+'\t'+str(len(nonshared_int))+'\t'+str(neo_int)
  with open(reportFile, "a") as f:
    f.write(out + '\n')
  
if __name__ == "__main__":
  parser = PDBParser(QUIET=True)
  io = PDBIO()
  
  code_id = sys.argv[1] #378132_Ag50-54_G107S.cdr.pdb
  native_path = sys.argv[2]  #'/content/gdrive/MyDrive/Ag_like/p53R175H/ranking/'+pdb_id+'.pdb'
  cdr_chain = sys.argv[3]
  cdr_list = sys.argv[4]
  cdr_resid = extractResid(cdr_list)
  ag_chain = sys.argv[5]
  ag_list = sys.argv[6]
  ag_resid = extractResid(ag_list)
  epi = sys.argv[7]
  epifrag_list = sys.argv[8]
  epifrag_resid = extractResid(epifrag_list)
  fullepitope_path = sys.argv[9]
  neo_resid = int(sys.argv[10])
  qdir = sys.argv[11]
  reportFile = qdir+'/ranking.txt'
  #wtepitope_path = sys.argv[12]
  
  native_target = dict()
  for native_resid, target_resid in zip(ag_resid, epifrag_resid):
    native_target[native_resid] = target_resid
  ag_contact_native = findSASAcontact(native_path, cdr_chain, cdr_resid, ag_chain, ag_resid)
  converted_ag_contact_native = [native_target[resid] for resid in ag_contact_native]
  print("Native",converted_ag_contact_native)
  
  opt_struct = parser.get_structure('opt', qdir+'/pdb/'+code_id+'.cdr.pdb')
  opt_cdr_resid = [res.id[1] for res in Selection.unfold_entities(opt_struct[0]["B"],"R") if res.id[1] in cdr_resid]
  target_struct = parser.get_structure('target', fullepitope_path)
  target_struct[0].add(opt_struct[0]["B"])
  io.set_structure(target_struct[0])
  io.save(qdir+'/pdb/'+code_id+'.EpiAb.pdb')
  ag_contact_target = findSASAcontact(qdir+'/pdb/'+code_id+'.EpiAb.pdb', "B", opt_cdr_resid, "A", [res.id[1] for res in target_struct[0]["A"].get_residues()])
  os.remove(qdir+'/pdb/'+code_id+'.EpiAb.pdb')
  ag_contact_epi = set(native_target.values())&set(ag_contact_target)
  print("Target",ag_contact_target)
  shared_int = set(converted_ag_contact_native)&set(ag_contact_target)
  print("Shared",shared_int, len(shared_int))
  nonshared_int = ag_contact_epi-set(converted_ag_contact_native)
  print("Nonshared",nonshared_int, len(nonshared_int))
  if neo_resid in ag_contact_target and neo_resid in converted_ag_contact_native:
    neo_int = True
    print("Found neo-site interaction!")
  else:
    neo_int = False
  out = code_id+'\t'+epi+'\t'+str(len(shared_int))+'\t'+str(len(nonshared_int))+'\t'+str(neo_int)
  with open(reportFile, "a") as f:
    f.write(out + '\n')