#!/usr/bin/env python

import numpy
from Bio.PDB import *
from Bio.PDB.SASA import ShrakeRupley
import sys
import os
import glob

def computeSASA(struct, computeRes):
  sr = ShrakeRupley()
  if computeRes == False:
    sr.compute(struct, level="S")
    return round(struct.sasa,2)
  sr.compute(struct, level="S")
  struct_sasa = round(struct.sasa,2)
  sr.compute(struct, level="R")
  sasa_res = dict()
  for res in computeRes:
    sasa_res[res.id[1]] = round(struct[0][res.get_parent().id][res.id].sasa,2)
  return struct_sasa, sasa_res

if __name__ == "__main__":
  parser = PDBParser(QUIET=True)
  io = PDBIO()
  sr = ShrakeRupley()

  qdir = sys.argv[1]
  opt_tcr_code = sys.argv[2]
  neo_resi = int(sys.argv[3])
  tcr_code = qdir.split("/")[-1]
  tcrdir = "./TCRgraft/"+"/".join(qdir.split("/")[2:-1])
  pdb = qdir.split("/")[3]
  tcr_chain = pdb.split("_")[-1]
  neoint = False
  
  #Before optimization
  tcr_struct = parser.get_structure('original', tcrdir+'/'+tcr_code+'.pdb')
  SASA_ogtcr = computeSASA(tcr_struct.copy(), computeRes=False)
  pmhc_struct = parser.get_structure('pmhc', tcrdir+'/target_aligned.pdb')
  pmhc_chain = list(pmhc_struct[0].get_chains())[0].id
  SASA_pmhc, SASA_neo = computeSASA(pmhc_struct.copy(), computeRes=[res for res in list(pmhc_struct[0][pmhc_chain].get_residues()) if res.id[1] in range(neo_resi-2,neo_resi+2)])
  complex_struct = pmhc_struct.copy()
  complex_struct[0].add(tcr_struct[0][tcr_chain])
  SASA_ogcom, SASA_neo_com_og = computeSASA(complex_struct.copy(), computeRes=[res for res in list(complex_struct[0][pmhc_chain].get_residues()) if res.id[1] in range(neo_resi-2,neo_resi+2)])
  BSA_ogcom = (SASA_pmhc+SASA_ogtcr)-SASA_ogcom
  #print("pMHC-SASA", round(SASA_pmhc,2), "\nogTCR-SASA", round(SASA_ogtcr,2), "\nComplex-SASA", round(SASA_ogcom,2), "\nBSA", round(BSA_ogcom,2))
  del complex_struct
  #print("Before deltaSASA-neosite", sum(SASA_neo.values())-sum(SASA_neo_com_og.values()))
  
  #After optimization
  opt_tcr = parser.get_structure('optimized', qdir+'/'+opt_tcr_code+'.pdb')
  SASA_tcr = computeSASA(opt_tcr.copy(), computeRes=False)
  complex_struct = pmhc_struct.copy()
  complex_struct[0].add(opt_tcr[0][tcr_chain])
  SASA_com, SASA_neo_com = computeSASA(complex_struct.copy(), computeRes=[res for res in list(complex_struct[0][pmhc_chain].get_residues()) if res.id[1] in range(neo_resi-2,neo_resi+2)])
  BSA_com = (SASA_pmhc+SASA_tcr)-SASA_com
  print(opt_tcr_code+"-SASA", "BSA", round(BSA_com,2))
  if sum(SASA_neo.values())-sum(SASA_neo_com.values()) > sum(SASA_neo.values())-sum(SASA_neo_com_og.values()):
    neoint = True
    print("New Candidate :-)", "\nImproved Neosite contact:", neoint,"(", "Before=", round(sum(SASA_neo_com_og.values()),2), "After=", round(sum(SASA_neo_com.values()),2),")")
  del complex_struct
  #print("After-complex", SASA_neo_com, sum(SASA_neo_com.values()))
  
  #Encounter testing
  encounter_results = dict()
  for enc_path in glob.glob("/".join(qdir.split("/")[:-1])+"/*.pdb"):
    encounter_results[enc_path.split("/")[-1].split(".")[0]] = True
    enc_struct = parser.get_structure('encounter', enc_path)
    SASA_enc = computeSASA(enc_struct.copy(), computeRes=False)
    complex_struct = enc_struct.copy()
    complex_struct[0].add(tcr_struct[0][tcr_chain])
    SASA_ogcom_enc = computeSASA(complex_struct.copy(), computeRes=False)
    BSA_ogcom_enc = (SASA_enc+SASA_ogtcr)-SASA_ogcom_enc
    del complex_struct
    complex_struct = enc_struct.copy()
    complex_struct[0].add(opt_tcr[0][tcr_chain])
    SASA_com_enc = computeSASA(complex_struct.copy(), computeRes=False)
    BSA_com_enc = (SASA_enc+SASA_tcr)-SASA_com_enc
    del complex_struct
    if BSA_com_enc < BSA_ogcom_enc:
      print("Decline", enc_path.split("/")[-1].split(".")[0])
      encounter_results[enc_path.split("/")[-1].split(".")[0]] = False

  if os.path.exists("/".join(qdir.split("/")[:-2])+"/optimizedTCRm_candidates.tsv"):
    with open("/".join(qdir.split("/")[:-2])+"/optimizedTCRm_candidates.tsv", "a") as File:
      File.write(opt_tcr_code+"\t"+pdb+"\t"+str(round(BSA_com,2))+"\t"+str(neoint)+"\t"+str("\t".join(str(v) for v in encounter_results.values()))+"\t"+qdir+'/'+opt_tcr_code+'.pdb\n')
  else:
    with open("/".join(qdir.split("/")[:-2])+"/optimizedTCRm_candidates.tsv", "w") as File:
      File.write("ID\tTCR\tBSA\tImporvedNeoContact\t"+str("\t".join(k for k in encounter_results.keys()))+"\tPath\n")
      File.write(opt_tcr_code+"\t"+pdb+"\t"+str(round(BSA_com,2))+"\t"+str(neoint)+"\t"+str("\t".join(str(v) for v in encounter_results.values()))+"\t"+qdir+'/'+opt_tcr_code+'.pdb\n')
