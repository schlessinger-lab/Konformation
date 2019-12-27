#!/usr/bin/python

import re,sys
import numpy as np
from x_helix_axis import *
from aa_residue import *
from Bio import SeqIO


# pocket   = [ [coordinates], sphere radius ]
pock_out   = [ [ 4.52,  7.98, 12.12], 3.3]  # 1
pock_hinge = [ [ 3.70,  9.76,  2.52], 2.7]  # 2
pock_sugar = [ [ 9.52,  8.59,  3.43], 2.5]  # 3
pock_up    = [ [12.20, 13.36,  0.14], 2.5]  # 4
pock_I5    = [ [ 6.97,  9.44, -2.34], 2.5]  # 5
pock_I5in  = [ [ 6.04,  8.94, -7.18], 2.5]  # 6
pock_I5dp  = [ [ 6.57, 12.34, -8.19], 2.5]  # 7
pock_dfg   = [ [ 6.28,  5.14, -6.99], 2.5]  # 8
pock_bp    = [ [ 8.88,  3.30,-10.83], 2.5]  # 9
Pockets = [ pock_out,  pock_hinge, pock_sugar, pock_up, pock_I5,
            pock_I5in, pock_I5dp,  pock_dfg,   pock_bp ]
No_Pock = [ None, None, None, None,
            [None,None,None,None,None,None,None,None,None] ]

##########################################################################
# Search thru PDB for HETATM ligands (exclude HOH, salts, additives, unnatAA),
# then calculate Ligand-Pocket distances if PDB has ligand(s)
def LigPocketOccupy( pdb_obj ):

  pdb_id  = pdb_obj.get_id()
  res_obj = pdb_obj.get_residues()

  LP_Dists = []
  for res in res_obj:
    Tag      = res.get_id() # ' ' = ATOM, 'H_' = HETATM, 'W' = water
    res_type = Tag[0]
    res_id   = Tag[1]
    if re.search(r'H_', res_type):
      res_name = res_type.split('_')[1]
      # Ignore hydrogens and salt/additive/unnaturalAA heteroatoms
      if not SaltAdditive(res_name) and not CheckUnnaturalAA(res_name):
        HA = [a for a in res if not re.search(r'H', a.get_name())]
        res_ha  = len(HA)	# number of heavy atom
        LP_Dists.append( [ pdb_id, res_name, res_id, res_ha, 
                           LigPockDistances(res) ] )

  if len(LP_Dists) == 0:
    LP_Dists.append(No_Pock)
  return LP_Dists


##########################################################################
# Calculate heavy-atom Ligand-Pocket distances, pocket coordinates are 
# pre-defined
def LigPockDistances( res ):

  Pocket_Occ = [None,None,None,None,None,None,None,None,None]

  for idx, pock in enumerate(Pockets):
    for atom in res:
      a_id = atom.get_name()
      if re.search(r'H', a_id):
        continue
      if Distance(res[a_id].coord, pock[0]) <= pock[1]:
        Pocket_Occ[idx] = True
        
  return Pocket_Occ


##########################################################################
# determine ligand type based on pharmacophere occupancy. Require ligand 
# presence
def LigandTypeDefinition(Poc, dfg_st):

# Pockets = [ pock_out (0),  pock_hinge (1), pock_sugar (2), pock_up (3),
#             pock_I5 (4), pock_I5in (5), pock_I5dp (6),  pock_dfg (7),
#             pock_bp (8) ]

  ligtype = None

  # for DFG-in...
  if dfg_st == 'in':
    # if no hinge fragment...
    if Poc[1] is None:
      # but with backpocket frag...
      if Poc[6] is True or Poc[7] is True or Poc[8] is True:
        lig = 'type-III'
      else:                 # nothing in the binding site
        lig = 'type-IV'
    # has hinge fragment...
    elif Poc[7] is True:    # but has dfg frag...
      if Poc[5] is True:
        lig= 'type-I1.5A'
      else:
        lig = 'other'
    elif Poc[8] is True:    # but has backpocket frag...
      lig = 'other'
    elif Poc[5] is True or Poc[6] is True:  # has I1.5 deep fragment...
      lig = 'type-I1.5A'
    elif Poc[4] is True and (Poc[5] is None or Poc[6] is None):   # no I1.5 frag
      lig = 'type-I1.5B'
    else:
      lig = 'type-I'

  # DFG-out and intermediates
  else:
    # If no hinge fragment...
    if Poc[1] is None:
      if Poc[7] is True and Poc[5] is True:    # and with dfg frag...
        lig = 'type-IIA'
      elif Poc[8] is True:  # but has dfg-pocket frag...
        lig = 'type-III'
      elif Poc[5] is True or Poc[6] is True:    # with I1.5 frag but no dfg...
        lig = 'type-I1.5A'  # few odd cases, 2XIR, 3C4C, lig points to dfg pocket
      else:
        lig = 'type-IV'     # nothing in the binding site
    # has hinge fragment
    else:
      if Poc[7] is True:    # has dfg frag...
        if Poc[8] is True:  # but with backpocket frag...
          lig = 'type-IIA'
        else:
          lig = 'type-IIA'
      elif Poc[8] is True:  # no dfg but backpocket frag...
        lig = 'type-III'
      elif Poc[5] is True or Poc[6] is True:    # no dfg but I1.5 frag...
        lig = 'type-I1.5A'	# odd cases, 2XIR, 3C4C, 4AG8, 5CSX
      else:                 # lig never reached deep
        lig = 'type-IIB'
        
  return lig


##########################################################################
def DescriptLigands(Ref, Tgt_List, Data):
# PDB_Coords=[pdb_name, H_Crds, N_Crds, C_Crds, G_Crds, R_Crds, T_Crds, L_Pock]
# L_Pock = [ [ pdb_id, lig_name, lig_id, lig_ha,
#              [pock_out, pock_hinge,pock_sugar,pock_up, pock_I5
#               pock_I5in,pock_IIb,  pock_dfg,  pock_bp           ]] ... ]

  Tgt_List.append(Ref)

  for Tgt in Tgt_List:

    pdb_id = Tgt[0]

    # pass on if no ligand (no Tgt[7][0][0]) is registered
    if Tgt[7][0][0] is None:
      Data[pdb_id]['lig_type'] = 'apo'
      continue

    # Prefer non-nucleotide ligands if 2 or more ligands are found.
    # only report the largest ligand if multiple ligands in 1 PDB
    # cannot distinguish if ligands are in binding site or not
    if len(Tgt[7]) > 1:
      Lig = [ lig for lig in Tgt[7] if not Nucleotides(lig[1]) ]
      if len(Lig) > 1:
        Lig.sort(key=lambda x: x[3], reverse=True)
        Data[pdb_id]['lig_multi'] = zip(*Lig)[1]
      Tt = Lig[0]
    else:
      Tt = Tgt[7][0]

    Data[pdb_id]['lig_name']   = Tt[1]
    Data[pdb_id]['lig_id']     = Tt[2]
    Data[pdb_id]['lig_ha']     = Tt[3]
    Data[pdb_id]['pock_out']   = Tt[4][0]
    Data[pdb_id]['pock_hinge'] = Tt[4][1]
    Data[pdb_id]['pock_sugar'] = Tt[4][2]
    Data[pdb_id]['pock_up']    = Tt[4][3]
    Data[pdb_id]['pock_I5']    = Tt[4][4]
    Data[pdb_id]['pock_I5in']  = Tt[4][5]
    Data[pdb_id]['pock_I5dp']  = Tt[4][6]
    Data[pdb_id]['pock_dfg']   = Tt[4][7]
    Data[pdb_id]['pock_bp']    = Tt[4][8]

    dfg_state = Data[pdb_id]['dfg_st']
    Data[pdb_id]['lig_type']   = LigandTypeDefinition(Tt[4], dfg_state)


##########################################################################
#
#   Peter M.U. Ung  @ MSSM
#
#   v1.0    17.08.18
#   v2.0    17.09.08 - change type-I.5B definition, check _II_other to _I.5A
#   v3.0    17.10.07   change type-III condition error
#   v4.0    18.03.15   add function to recognize multiple ligands and select
#                      most important one (largest non-nucleotide) 
# 
##########################################################################
