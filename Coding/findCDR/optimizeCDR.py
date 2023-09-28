#!/usr/bin/env python

from Bio.PDB import *
import numpy
import os
import sys
import glob
from multiprocessing import Process

def getClosestRes(tar_res, match_atom):
  tar_ca = tar_res["CA"]
  tar_distance = [tar_ca - match_ca for match_ca in match_atom]
  min_tar_distance_ind = tar_distance.index(min(tar_distance))
  return match_atom[min_tar_distance_ind].get_parent()

def findBetterRes(match_struct, og_struct, ab_opt):
  match_atom = [atom for atom in match_struct[0].get_atoms() if atom.name=="CA"]
  for ag_res in og_struct[0]["A"].get_residues():
    ag_inMatch = getClosestRes(og_struct[0]["A"][ag_res.id[1]], match_atom)
    print(ag_res.id[1], ag_res.resname, ag_inMatch.resname)
    if ag_res.resname != ag_inMatch.resname:
      return numpy.nan, numpy.nan, numpy.nan
  ab_inMatch = getClosestRes(og_struct[0]["B"][ab_opt], match_atom)
  return ab_inMatch.id[1], ab_inMatch.resname, ab_inMatch.get_full_id()[2]

def checkDoAc(atom, type):
  if type == "Donor":
    H_dict = {"All":["N"], "SER":["OG"], "THR":["OG1"], "LYS":["NZ"], "ASN":["ND2"], "GLN":["NG2"]}
  elif type == "Acceptor":
    H_dict = {"All":["O","OXT"], "ASP":["OD1","OD2"], "ASN":["OD1"], "GLU":["OG1","OG2"], "GLN":["OG1"]}
  allow = 0
  if atom.get_parent().resname in H_dict.keys():
    if atom.name in H_dict[atom.get_parent().resname]:
      allow += 1
  else:
    if atom.name in H_dict["All"]:
      allow += 1
  return allow

def findClashContact(Ab_res, ns, _cutoff_nsdist, _cutoff_overlap, _cutoff_contact, _Hbond_allowance):
  dist_dict = {"C":1.700,"N":1.625,"O":1.480,"S":1.782,"H":1.000}
  result = {"Clash":[],"Allcontact":[]}
  for target in Ab_res.get_atoms(): 
    target_donor = checkDoAc(target, "Donor")
    target_acceptor = checkDoAc(target, "Acceptor")
    close_atoms = ns.search(target.coord, _cutoff_nsdist)
    for n, close_atom in enumerate(close_atoms):
      allowance=False
      close_donor = checkDoAc(close_atom, "Donor")
      close_acceptor = checkDoAc(close_atom, "Acceptor")
      if target_donor+close_acceptor == 2:
        allowance = True
      elif target_acceptor+close_donor == 2:
        allowance = True
      distance = target - close_atom
      for symbol in dist_dict.keys():
        if symbol in target.name:
          target_radii = dist_dict[symbol]
        if symbol in close_atom.name:
          close_atom_radii = dist_dict[symbol]
      overlap_dist = (target_radii+close_atom_radii)-distance
      if allowance == True:
        if overlap_dist-_Hbond_allowance < _cutoff_overlap:
          if overlap_dist >= _cutoff_contact:
            if close_atom.get_parent().id[1] not in result["Allcontact"]:
              result["Allcontact"].append(close_atom.get_parent().id[1])
        else:
          if close_atom.get_parent().id[1] not in result["Clash"]:
            result["Clash"].append(close_atom.get_parent().id[1])
      else:
        if overlap_dist < _cutoff_overlap:
          if overlap_dist >= _cutoff_contact:
            if close_atom.get_parent().id[1] not in result["Allcontact"]:
              result["Allcontact"].append(close_atom.get_parent().id[1])
        else:
          if close_atom.get_parent().id[1] not in result["Clash"]:
            result["Clash"].append(close_atom.get_parent().id[1])
  return result

def saveOptimizedCDR(match_path, og_path, ns, _cutoff_nsdist, _cutoff_overlap, _cutoff_contact, _Hbond_allowance, qdir):
  code_standard = {
  'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
  'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
  'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
  'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G'}
  cdr_code = og_path.split("contact")[0].split("_EpiAb_")[0].split("/")[-1]
  
  og_struct = parser.get_structure('original', og_path)
  ab_opt = int(og_path.split("contact")[0].split("_")[-1])
  ab_res = og_struct[0]["B"][ab_opt].resname
  
  print(match_path)
  match_struct = parser.get_structure('match', match_path)
  betterRes_id, betterRes, betterRes_chain = findBetterRes(match_struct, og_struct, ab_opt)
  if numpy.isnan(betterRes_id) == False:
    print(betterRes_chain, betterRes, betterRes_id)
    cdr_res = match_struct[0][betterRes_chain][betterRes_id]
    result = findClashContact(cdr_res, ns, _cutoff_nsdist, _cutoff_overlap, _cutoff_contact, _Hbond_allowance)
    if len(result['Clash']) == 0:
      if betterRes == ab_res:
        io.set_structure(og_struct[0]["B"])
        io.save(qdir+'/'+cdr_code+'_'+str(ab_opt)+'_optimized.pdb')
      elif betterRes != "PRO" or betterRes != "CYS":
        new_og_struct = og_struct.copy()
        cdr_res.id = tuple((' ', 10000000, ' '))
        insert_pos = [res for res in new_og_struct[0]["B"]].index(new_og_struct[0]["B"][ab_opt])
        new_og_struct[0]["B"].detach_child((' ', ab_opt, ' '))
        new_og_struct[0]["B"].insert(insert_pos, cdr_res)
        new_og_struct[0]["B"][10000000].id = tuple((' ', ab_opt, ' '))
        io.set_structure(new_og_struct[0]["B"])
        io.save(qdir+'/'+cdr_code+'_'+code_standard[ab_res]+str(ab_opt)+code_standard[betterRes]+'_optimized.pdb')
      print("Optimized-CDR saved")
  os.remove(match_path)

if __name__ == "__main__":
  parser = PDBParser(QUIET=True)
  io = PDBIO()
  _cutoff_nsdist = 6
  _cutoff_overlap = 0.6
  _Hbond_allowance = 0.4
  _cutoff_contact = -0.4
  
  qdir = sys.argv[1]
  Gpath = sys.argv[2]
  target_path = sys.argv[3]
  
  target = parser.get_structure('fullepitope', target_path)
  atom_list = [_ for _ in target[0]["A"].get_atoms()]
  ns = NeighborSearch(atom_list)
  
  og_path = glob.glob(qdir+'/*contact.pdb')[0]
  for match_path in glob.glob(qdir+'/'+Gpath+'.match.struct/match**.pdb'):
    process = Process(target=saveOptimizedCDR, args=(match_path, og_path, ns, _cutoff_nsdist, _cutoff_overlap, _cutoff_contact, _Hbond_allowance, qdir))
    process.start()
    process.join()
