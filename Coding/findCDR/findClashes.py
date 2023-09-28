#!/usr/bin/env python

from Bio.PDB import *
import numpy
import os
import sys

def extractPos(subline):
  chain = subline.split("[")[1].split("'")[1]
  pos = list()
  for ele in subline.split("[")[2:]:
    res=list()
    try:
      for subele in ele.split(', '):
        res.append((int(subele)))
    except:
      res.append((int(subele.split(']')[0])))
    if len(res) > 1:
      pos.extend([n for n in range(int(res[0]),int(res[1])+1)])
    else:
      pos.extend([res[0]])
  return chain, pos
  
class ResSelect(Select):
    def __init__(self,select_res,chain_id):
        self.select_res = select_res
        self.chain_id = chain_id
    def accept_residue(self, res):
        if res.id[1] in self.select_res and res.parent.id == self.chain_id:
            return 1
        else:
            return 0

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = numpy.zeros((len(chain_one), len(chain_two)))
    for row, residue_one in enumerate(chain_one):
        for col, residue_two in enumerate(chain_two):
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer

def extractAgPDB(matchFile):
  match_struc = parser.get_structure('match', matchFile)
  match_chains = list(match_struc[0].get_chains())
  Ag_pos = list(res.id[1] for res in match_struc[0][match_chains[0].id].get_residues())
  #Ag_pos = sorted([res.id[1] for res in match_struc[0][match_chains[0].id]])
  return match_chains[0].id, Ag_pos

def findClash(ind, pdb, Ab_chain, Ab_pos, Ag_chain, Ag_pos, epi_pos, fullEpiFile, saveDir):
    full_struc = parser.get_structure('full', "/gs/hs0/tga-sca/FAH/all_PDB/ncAA/"+pdb+".pdb")
    io.set_structure(full_struc[0])
    io.save(saveDir+"/"+ind+"_Ag"+str(Ag_pos[0])+"-"+str(Ag_pos[-1])+'_Ab.pdb', ResSelect(Ab_pos,Ab_chain))
    Ab_struc = parser.get_structure('Ab', saveDir+"/"+ind+"_Ag"+str(Ag_pos[0])+"-"+str(Ag_pos[-1])+'_Ab.pdb')
    Ab_chains = list(Ab_struc[0].get_chains())
    Ab_chains[0].id = 'B'

    io.set_structure(full_struc[0])
    io.save(saveDir+"/"+ind+"_Ag"+str(Ag_pos[0])+"-"+str(Ag_pos[-1])+'_Ag.pdb', ResSelect(Ag_pos,Ag_chain))
    Ag_struc = parser.get_structure('Ag', saveDir+"/"+ind+"_Ag"+str(Ag_pos[0])+"-"+str(Ag_pos[-1])+'_Ag.pdb')
    Ag_chains = list(Ag_struc[0].get_chains())
    Ag_chains[0].id = 'G'
    
    Ag_struc[0].add(Ab_chains[0])
    alt_atoms = []
    for alt_res in Ag_struc[0]['G']:
      if alt_res.id[1] in Ag_pos:
        alt_atoms.append(alt_res['CA'])    
    epitope_struc = parser.get_structure('epitope', fullEpiFile)
    ref_atoms = []
    for ref_res in epitope_struc[0]['A']:
      if ref_res.id[1] in epi_pos:
        ref_atoms.append(ref_res['CA'])
    try:
      sup.set_atoms(ref_atoms, alt_atoms)
      sup.apply(Ag_struc[0].get_atoms())
    except:
      print("Alignment Error:", ind, pdb, Ab_chain, Ab_pos, Ag_chain, Ag_pos, epi_pos)
      return True
    io.set_structure(Ag_struc[0])
    io.save(saveDir+"/"+ind+"_Ag"+str(Ag_pos[0])+"-"+str(Ag_pos[-1])+'_AbAg_aligned.pdb')
    
    AbAg_align = parser.get_structure('align', saveDir+"/"+ind+"_Ag"+str(Ag_pos[0])+"-"+str(Ag_pos[-1])+'_AbAg_aligned.pdb')
    dis_matrix = calc_dist_matrix(AbAg_align[0]['B'], epitope_struc[0]['A'])
    print(ind, pdb, Ab_chain, Ab_pos, Ag_chain, Ag_pos, epi_pos)
    print(numpy.min(dis_matrix))
    if round(numpy.min(dis_matrix),1) < 4.0: #distance between alpha-carbon atom pair of the connected residues ~ 3.8
      print("Steric clash!")
      os.remove(saveDir+"/"+ind+"_Ag"+str(Ag_pos[0])+"-"+str(Ag_pos[-1])+'_AbAg_aligned.pdb')
      return True
    else:
      print("OK!")
      savename = saveDir+"/"+ind+"_Ag"+str(Ag_pos[0])+"-"+str(Ag_pos[-1])+"_EpiAb"
      if os.path.exists(savename) == False:
        os.mkdir(savename)
        epitope_struc[0].add(Ag_struc[0]['B'])
        io.set_structure(epitope_struc[0])
        io.save(savename+"/"+ind+"_Ag"+str(Ag_pos[0])+"-"+str(Ag_pos[-1])+'_EpiAb.pdb')
      os.remove(saveDir+"/"+ind+"_Ag"+str(Ag_pos[0])+"-"+str(Ag_pos[-1])+'_AbAg_aligned.pdb')
      return False

if __name__ == "__main__":
  parser = PDBParser(QUIET=True)
  io = PDBIO()
  sup = Superimposer()
  
  epifrag = sys.argv[1]
  matchFile = sys.argv[2]
  AbAgfile = sys.argv[3]
  fullEpiFile = sys.argv[4]
  saveDir = sys.argv[5]
  reportFile = saveDir+"/CDR-like_candidates.lst"
  
  epi_pos1 = int(epifrag.split("_")[2].split("-")[0])
  epi_pos2 = int(epifrag.split("_")[2].split("-")[1])
  epi_pos = [n for n in range(epi_pos1,epi_pos2+1)]
  
  Ag_chain, Ag_pos = extractAgPDB(matchFile)
  with open(AbAgfile, "r") as file:
    for line in file:
      ind = str(line.split()[0])
      pdb = str(line.split()[2])
      Ab_chain, Ab_pos = extractPos(line.split("]]")[0].split("\t")[-1])
      if findClash(ind, pdb, Ab_chain, Ab_pos, Ag_chain, Ag_pos, epi_pos, fullEpiFile, saveDir) == False:
        with open(reportFile, "a") as outfile:
          out=str(ind+'\t'+pdb+'\t'+Ab_chain+'\t'+str(Ab_pos)+'\t'+Ag_chain+'\t'+str(Ag_pos)+'\t'+epifrag+'\t'+str(epi_pos))
          outfile.write(out+'\n')
      os.remove(saveDir+"/"+ind+"_Ag"+str(Ag_pos[0])+"-"+str(Ag_pos[-1])+'_Ab.pdb')
      os.remove(saveDir+"/"+ind+"_Ag"+str(Ag_pos[0])+"-"+str(Ag_pos[-1])+'_Ag.pdb')


    
