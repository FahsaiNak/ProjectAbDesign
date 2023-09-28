#!/usr/bin/env python

from Bio.PDB import *
from biopandas.pdb import PandasPdb
import sys
from Bio import SeqIO
import pandas as pd

code_standard = {
    'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
    'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
    'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
    'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G'}
class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0
class ResSelect(Select):
    def __init__(self,select_res):
        self.select_res = select_res
    def accept_residue(self, res):
        if res.id[1] in self.select_res:
            return 1
        else:
            return 0
def insertResidues(tempchain, addchain):
  startpoint = list(tempchain.get_residues())[-1].id[1]
  for nres, res in enumerate(addchain.get_residues()):
    res.id = tuple((' ', startpoint+nres+1, ' '))
    tempchain.add(res)
  return tempchain

if __name__ == "__main__":
  parser = PDBParser(QUIET=True)
  io = PDBIO()
  sup = Superimposer()
  
  target_path = sys.argv[1]
  wk_dir = sys.argv[2] #"/gs/hs0/tga-sca/FAH/rosettaMHC/p53R175H/Ag-like/linear/neoCDRcandidates/shortlist/TCRgraft/TCRb3/${TCRpdb}"
  template_name = wk_dir.split("/")[-1]
  template_pdb = wk_dir.split("/")[-1].split("_")[0]
  template_b_chain = wk_dir.split("/")[-1].split("_")[1]
  
  #remove heteroatoms from a original TCR structure
  TCR_struct_og = parser.get_structure('TCRog', wk_dir+"/"+template_pdb+".pdb")
  io.set_structure(TCR_struct_og)
  io.save(wk_dir+"/"+template_pdb+"_cleaned.pdb", NonHetSelect())
  
  #extract the single variable domain from a TCR beta chain
  template_struct = parser.get_structure('fullTCR', wk_dir+"/"+template_pdb+"_cleaned.pdb")
  template_b_chain_seq = []
  for res in template_struct[0][template_b_chain]:
    template_b_chain_seq.append(code_standard[res.resname])
  template_b_chain_seq = ''.join(template_b_chain_seq)
  fasta_sequences = list(SeqIO.parse(open(wk_dir+"/"+template_name+"_VB.fasta"),'fasta'))
  template_vb_len = template_b_chain_seq.find(str(fasta_sequences[0].seq))+len(fasta_sequences[0].seq)
  io.set_structure(template_struct[0][template_b_chain])
  io.save(wk_dir+"/"+template_name+"_TCRvb.pdb", ResSelect([n for n in range(list(template_struct[0][template_b_chain].get_residues())[:template_vb_len][-1].id[1]+1)]))

  #align the target pMHC to the TCR-template pMHC orientation
  TCR_df = pd.read_csv('/gs/hs0/tga-sca/FAH/human_tcr_STCRDab/summary.tsv', sep='\t')
  tcrmhc_chain = TCR_df["mhc_chain1"][(TCR_df["pdb"] == template_pdb) & (TCR_df["Bchain"] == template_b_chain)].values[0]
  template_struct = parser.get_structure('fullTCR', wk_dir+"/"+template_pdb+"_cleaned.pdb")
  ref_atoms = []
  for ref_res in list(template_struct[0][tcrmhc_chain].get_residues())[:180]:
    ref_atoms.append(ref_res['CA'])  
  target_struct = parser.get_structure('full', target_path)
  alt_atoms = []
  for alt_res in list(target_struct[0][list(target_struct[0].get_chains())[0].id].get_residues())[:180]:
    alt_atoms.append(alt_res['CA'])
  sup.set_atoms(ref_atoms, alt_atoms)
  sup.apply(target_struct[0].get_atoms())
  io.set_structure(target_struct[0])
  io.save(wk_dir+"/target_aligned.pdb")
  
  TCRvb_struct = parser.get_structure('TCRvb', wk_dir+"/"+template_name+"_TCRvb.pdb")
  target_struct[0].add(TCRvb_struct[0][template_b_chain])
  io.set_structure(target_struct[0])
  io.save(wk_dir+"/target_pMHCTCRvb.pdb")
  
  #save MHC and pMHC structures from TCR template
  io.set_structure(template_struct[0][tcrmhc_chain])
  io.save(wk_dir+"/"+template_pdb+"_MHC.pdb", ResSelect([n for n in range(list(template_struct[0][tcrmhc_chain].get_residues())[:180][-1].id[1]+1)]))
  tcrag_chain = TCR_df["antigen_chain"][(TCR_df["pdb"] == template_pdb) & (TCR_df["Bchain"] == template_b_chain)].values[0]
  mhc_struct = parser.get_structure('mhc', wk_dir+"/"+template_pdb+"_MHC.pdb")
  for nres, res in enumerate(mhc_struct[0][tcrmhc_chain].get_residues()):
    res.id = tuple((' ', nres+1, ' '))
  pmhc_struct = insertResidues(mhc_struct[0][tcrmhc_chain], template_struct[0][tcrag_chain])
  io.set_structure(pmhc_struct)
  io.save(wk_dir+"/"+template_pdb+"_pMHC.pdb")
  
  TCRvb_struct = parser.get_structure('TCRvb', wk_dir+"/"+template_name+"_TCRvb.pdb")
  pmhc_struct = parser.get_structure('pMHC', wk_dir+"/"+template_pdb+"_pMHC.pdb")
  pmhc_struct[0].add(TCRvb_struct[0][template_b_chain])
  io.set_structure(pmhc_struct[0])
  io.save(wk_dir+"/"+template_pdb+"_pMHCTCRvb.pdb")
  
  