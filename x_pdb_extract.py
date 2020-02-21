#!/usr/bin/env python3

import re,os,sys
import numpy as np
from aa_residue import AA
#from x_helix_axis  import *
#from x_dfg_torsion import *
from x_fasta_parse import CheckSequence
from x_ligand_type import LigPocketOccupy

from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
p = PDBParser(PERMISSIVE=1)

##########################################################################
class ParsePDB(object):

  def __init__( self, h_seq=None, n_seq=None, c_seq=None,
                      g_seq=None, f_seq=None, t_seq=None,
                      pdb_dir=None, corr={}   ):
    self.h_seq   = h_seq    # Helix sequence, [1]
    self.n_seq   = n_seq    # N-dom sequence, [1]
    self.c_seq   = c_seq    # C-dom sequence, [1]
    self.f_seq   = f_seq    # DFG-F sequence, [1]
    self.g_seq   = g_seq    # Gate sequence, [1]

#    self.r_seq   = r_seq    # R-spine sequences, [4]
    self.t_seq   = t_seq    # C-spine sequences, [8]
    self.pdb_dir = pdb_dir  # default PDB directory, unless it is changed
    self.corr    = corr     # sequence correction list

  def __call__( self, pdb ):
    return self.extract_pdb( pdb )

########################################################
  def extract_pdb( self, pdb ):
    
    pdb_name = pdb.split('/')[-1]
    pdb_id   = pdb.split('/')[-1].split('.')[0]
    with open('_TEMP.missing.'+pdb_name, 'w') as missing:
      print('>>> Current PDB: '+pdb_id)

      # Skip if protein is not found in fasta library or sequence has blank 
      # residue '-' in the sequence
      if pdb_id not in self.h_seq:
        print('\n  \033[31m#2# PDB Skip:\033[0m Cannot find in FASTA library: '+pdb_id)
        return None
      if CheckSequence( self.h_seq[pdb_id] ) is False:
        print('\n  \033[31m#2# PDB Warning:\033[0m FASTA has missing residue: {0} - {1}'.format(
              pdb_id, str(self.h_seq[pdb_id]) ))
        return None

      # Found PDB in pdb_dir 
      for p_dir in self.pdb_dir.split(','):
        if not p_dir:
          continue
        if re.search(r'~', p_dir):
          sys.exit('\n  \033[31m#2# PDB FATAL:\033[0m Python does not recognize "~" home directory\n            Use full directory name: '+p_dir)
        if not os.path.isfile( pdb ):
          sys.exit('\n  \033[31m#2# PDB FATAL:\033[0m PDB not found in directory: '+pdb_name)

      pdb_obj = p.get_structure(pdb_name, pdb)

      L_Pock = LigPocketOccupy( pdb_obj )

      H_Crds = ExtractPDBCoords( pdb_obj, self.h_seq[pdb_id], 0 ) # (XHELIX)
      N_Crds = ExtractPDBCoords( pdb_obj, self.n_seq[pdb_id], 4 ) # (ZNDOM)
      C_Crds = ExtractPDBCoords( pdb_obj, self.c_seq[pdb_id], 4 ) # (ZCDOM)
      F_Crds = ExtractPDBCoords( pdb_obj, self.f_seq[pdb_id], 4 ) # (DFG_F)
      G_Crds = ExtractPDBCoords( pdb_obj, self.g_seq[pdb_id], 3 ) # (Gate)

      ## R- and C-spines residues are already coming in as list(dict of seq)
      R_Crds = [] # [ExtractPDBCoords(pdb_obj, seq[pdb_id]) for seq in self.r_seq]
      T_Crds = [] # [ExtractPDBCoords(pdb_obj, seq[pdb_id]) for seq in self.t_seq]
      
      # If the coordinates collection failed in the previous step, check
      # if correction data for failed residues is available for replacement,
      # otherwise output as None and ignore this PDB in future calculations
      # and marked as missing in an output file
      if H_Crds is None:
        if pdb_id in self.corr and self.corr[pdb_id][0] is not None:
          H_Crds = self.corr[pdb_id][0]
          print('  #1# PDB Info: Accepted coordinates correction: '+pdb_id+' Helix')
        else:
          missing.write(pdb_name+'|Helix|'+''.join(self.h_seq[pdb_id])+'\n')
#          return None
      if N_Crds is None:
        if pdb_id in self.corr and self.corr[pdb_id][1] is not None:
          N_Crds = self.corr[pdb_id][1]
          print('  #1# PDB Info: Accepted coordinates correction: '+pdb_id+' N_dom')
        else:
          missing.write(pdb_name+'|N_dom|'+''.join(self.n_seq[pdb_id])+'\n')
#          return None
      if C_Crds is None:
        if pdb_id in self.corr and self.corr[pdb_id][2] is not None:
          C_Crds = self.corr[pdb_id][2]
          print('  #1# PDB Info: Accepted coordinates correction: '+pdb_id+' C_dom')
        else:
          missing.write(pdb_name+'|C_dom|'+''.join(self.c_seq[pdb_id])+'\n')
#          return None
      if F_Crds is None:
        if pdb_id in self.corr and self.corr[pdb_id][3] is not None:
          F_Crds = self.corr[pdb_id][3]
          print('  #1# PDB Info: Accepted coordinates correction: '+pdb_id+' DFG_F')
        else:
          missing.write(pdb_name+'|DFG_F|'+''.join(self.f_seq[pdb_id])+'\n')
#          return None
      if G_Crds is None:
        if pdb_id in self.corr and self.corr[pdb_id][4] is not None:
          G_Crds = self.corr[pdb_id][4]
          print('  #1# PDB Info: Accepted coordinates correction: '+pdb_id+' Gate')
        else:
          missing.write(pdb_name+'|Gate|'+''.join(self.g_seq[pdb_id])+'\n')
#          return None


      return [pdb_id, H_Crds, N_Crds, C_Crds, G_Crds, F_Crds, T_Crds, L_Pock]


##########################################################################
## From the input Biopython PDB Object, extract the coordinates and info of
## 'res_num' rsidues corresponding to the supplied 'Query Sequence'
def ExtractPDBCoords( PDB, Query_Seqs, posit_backup ):

  # Quit searching if None is provided
  if Query_Seqs is None:
    print('\n  \033[31m#2# PDB Skip:\033[0m No input data for coordinate extraction: \033[31m{}\033[0m'.format(PDB.get_id()))
    return None

  for idx, Query_Seq in enumerate(Query_Seqs):
    res_num = len(Query_Seq)
    pdb_id  = PDB.get_id()
    Res_Obj = PDB.get_residues()
    if re.search(r'.pdb', pdb_id):
      pdb_id = pdb_id.split('.pdb')[0]
    print('>> Query Sequence:\t{0} -\t{1}'.format(pdb_id, ''.join(Query_Seq)))

    ## Convert BioPython Residue Object into List of Residues
    Residues = []
    for res in Res_Obj:
      if re.search(r"H_|W", res.get_id()[0]): continue

      resname   = res.get_resname()
      resid     = res.get_id()[1]
      bb_crds, ca_crd, cg_crd, avg_crd, cb_crd, cd_crd = ResidueCoords( res )

      Residues.append( [ resname, resid, bb_crds, ca_crd, cg_crd, 
                                          avg_crd, cb_crd, cd_crd ] )

    # Convert the target sequence into 3-letter AA name. Number of residue
    # depends on the input sequence length (variable)
    Target_Seq = [AA(Query_Seq[i]) for i in range(res_num)]

    Found = LocateTargetSeq( Target_Seq, Residues, pdb_id )
    if Found is not None:
      if idx == 0:
        return Found
#        break
      if idx == 1:
        print('  @ Found residue with 2nd Sequence: \033[31m{}\033[0m'.format(AA(Found[posit_backup][0])))
        return [ Found[posit_backup] ]
    else:
      if idx == 0:
        print('  \033[31m!! Warning:\033[0m Fail to find residue with Primary Sequence')
        print('  \033[31m!!\033[0m          Use Secondary recognition sequence')
        continue
      else:
        print('\n  \033[31m#2# PDB Skip:\033[0m Fail to find residue with both 1st/2nd Sequences: '+pdb_id)
        return None


##########################################################################

def LocateTargetSeq( Target_Seq, Residues, pdb_id ):

  # At least 5 residues are needed for a positive match
  while ( len(Target_Seq) >= 5 ):
    res_num = len(Target_Seq)
    # Iterate through the entire protein sequence and match the sequence to
    # the entire target sequence AA in 'Target_Seq' (vary in length)
    Found, matched = [], False
    for idx, residue in enumerate(Residues):
      # If the current position is reaching the end of the sequence and
      # only 'res_num' residue away, meaning won't have the exact matching
      # to the 'Target_Seq', break the search
      if (idx + res_num-1) >= len(Residues):
        break

      Found = []
      # If found a match resname, check the next residue for next in 
      # 'Target_Seq' until matching all resnames in 'Target_Seq'
      for step in range(res_num):
        if Residues[idx+step][0] == Target_Seq[step]:
          Found.append( Residues[idx+step] )
#          print(Residues[idx+step][0]+str(Residues[idx+step][1]))
        else:
          break   # If no match, break search and step to next residue

      if len(Found) == res_num:
        matched = True
        break
      else:
        continue

    # If no matching sequence is found after moving thru the entire protein,
    # run it one more time with truncation of 1 residue on both ends
    # if failed after rerun with 5 total residues, return None
    #
    # [ resname, resid, bb_crds, ca_crd, cg_crd, avg_crd, cb_crd ]
    if matched is False:
      Trunc_Seq  = Target_Seq[1:-1]
      Target_Seq = Trunc_Seq
      if len(Target_Seq) <= 3:
        print('\n  \033[31m#2# PDB Skip:\033[0m Cannot find sequence in\t{0}'.format(pdb_id))
        return None
      else:
        print('  \-33[31m## PDB:\033[0m Cannot find match in {0}. Shortened to {1}'.format(
                pdb_id, len(Target_Seq) ) )
    else:
      print(' Matched sequence in\t{0} -\033[31m\t{1}-{3}-{2}\033[0m'.format(
              pdb_id, Found[0][1], Found[-1][1],
              ''.join([AA(Found[i][0]) for i in range(0,len(Found))])  ) )
      return Found


##########################################################################
## Optional file contain correction for PDBs that have missing residues in
## the sequence matching region. Run thru the program first to identify
## which PDB has what missing. The correction residues should be in N-CA-C
## order and center on the key residues only
##   Format of corrected file: correct.<PDB file> -->correct.1ATP_E.1atp.pdb
def CoordCorrect( missing, pdb_dir ):

  dic = {}
  if missing is None: return dic

  # Read in correction file and load correction data into database
  with open(missing, 'r') as fi:
    for l in fi:
      pdb_name, typ, seq = l.split('|')
      pdb_id = pdb_name.split('.')[0]

      for p_dir in pdb_dir.split(','):
        if not p_dir:
          continue
        if re.search(r'~', p_dir):
          sys.exit('\n  \033[31m#2# PDB FATAL:\033[0m Python does not recognize "~" home directory\n            Use full directory name: '+p_dir)
        if os.path.isfile(p_dir+'/correct.'+typ+'.'+pdb_name):
          pdb = p_dir+'/correct.'+typ+'.'+pdb_name
      if pdb is None:
        print('# No correction for {0} : {1}'.format(pdb_id, typ))
        continue
      px      = p.get_structure(pdb_id, pdb)
      All_Res = px.get_residues()

      Residue = []
      for res in All_Res:
        resname = res.get_resname()
        resid   = res.get_id()[1]
        bb_crds, ca_crd, cg_crd, avg_crd, cb_crd, cd_crd = ResidueCoords( res )
        Residue.append( [ resname, resid, bb_crds, ca_crd, cg_crd, 
                                          avg_crd, cb_crd ] )

      # Put correction data into separated arrays,
      # [Helix, N-Dom, C-Dom, DFG_F, Gate, Rs1, Rs2, Rs3, Rs4, 
      #  Cs1, Cs2, Cs3, Cs4, Cs5, Cs6, Cs7, ]
      if pdb_id not in dic:
        dic[pdb_id] = [ None, None, None, None, None, None, None, None,
                        None, None, None, None, None, None, None, ]

      if   typ == 'Helix': dic[pdb_id][0]  = Residue
      elif typ == 'N_dom': dic[pdb_id][1]  = Residue
      elif typ == 'C_dom': dic[pdb_id][2]  = Residue
      elif typ == 'DFG_F': dic[pdb_id][3]  = Residue
      elif typ == 'Gate':  dic[pdb_id][4]  = Residue
#      elif typ == 'Rs1':   dic[pdb_id][5]  = Residue
#      elif typ == 'Rs2':   dic[pdb_id][6]  = Residue
#      elif typ == 'Rs3':   dic[pdb_id][7]  = Residue
#      elif typ == 'Rs4':   dic[pdb_id][8]  = Residue
#      elif typ == 'Cs1':   dic[pdb_id][9]  = Residue
#      elif typ == 'Cs2':   dic[pdb_id][10] = Residue
#      elif typ == 'Cs3':   dic[pdb_id][11] = Residue
#      elif typ == 'Cs4':   dic[pdb_id][12] = Residue
#      elif typ == 'Cs5':   dic[pdb_id][13] = Residue
#      elif typ == 'Cs6':   dic[pdb_id][14] = Residue
#      elif typ == 'Cs7':   dic[pdb_id][15] = Residue
#      elif typ == 'Cs8':   dic[pdb_id][16] = Residue

  return dic


##########################################################################
## Extract the coordinates of the residue -- backbone N,CA,C, CA, and average
## of sidechain from CB
def ResidueCoords( res ):

  ca_coord  = None
  cb_coord  = None
  cd_coord  = None
  cg_coord  = None
  bb_coords = []
  sc_coords = []

  for atom in res:
    a_id = atom.get_name()
    if re.search(r'H', a_id):   # Ignore hydrogen; fail: Arg NH1 NH2 nitrogen
      continue
    else:
      # Coordinates of backbone atoms N,CA,C
      if a_id == 'CA' or a_id == 'C' or a_id == 'N' or a_id == 'H':
        bb_coords.append(res[a_id].get_coord())
      # Coordinates of sidechain atoms
      else:
        sc_coords.append(res[a_id].get_coord())

      # Coordinates of other atoms: CA, CB, CG
      # Use CA as representation of mainchain of AA
      # Use CG as representation of side chain; most AA has atoms beyond CG
      # Retreive CB for DFG-vector calculation and cg_vec check
      if res.has_id('CA'):
        ca_coord = res['CA'].get_coord()
      if res.has_id('CB'):
        cb_coord = res['CB'].get_coord()
      if res.has_id('CD'):
        cd_coord = res['CD'].get_coord()

      if res.has_id('CG'):
        cg_coord = res['CG'].get_coord()
      elif res.has_id('CG2'): # Ile, Thr, Val has branch at CB
        cg_coord = res['CG2'].get_coord()

  avg_coord = np.mean(sc_coords, axis=0)

  return bb_coords, ca_coord, cg_coord, avg_coord, cb_coord, cd_coord 


##########################################################################
#
#   Peter M.U. Ung  @ MSSM
#   
#   v0.1    17.01.30
#   v0.2    17.02.01
#   
#   v1.0    17.02.02
#   v2.0    17.03.10
#   v3.0    17.03.19    check multiple pdb_directory
#   v4.0    17.04.12    bug fix; off-load the NoneType handling to subprocesses
#                       instead of ignoring all of them
#   v5.0    17.05.08    differentiate CG and CB coordinates/vectors
#   v6.0    17.08.18    include ligands during extraction
#   v6.5    17.10.09    if CB is missing, use CA as substitute
#   v7.0    18.03.13    updated the codes to accept PDB with full path
#   v8.0    18.03.15    masked all R/C-spine extraction, add 2nd recognition
#                       sequence search function
#
##########################################################################
