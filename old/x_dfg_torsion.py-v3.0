#!/usr/bin/python


import re,os,glob,sys
import numpy as np
import multiprocessing
from aa_residue    import *
from x_helix_axis  import *
from CommonUtility import *
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
p = PDBParser(PERMISSIVE=1)
np.seterr(invalid='ignore')

##########################################################################
##
def DFGTorsionAngle( Ref_Coords, Tgt_Coords, Data, output ):

  # Input_Coords = [pdb_name, H_Crds, N_Crds, C_Crds, G_Crds, R_Crds, T_Crds]
  #     x_Coords = [resname, resid, bb_crds, ca_crd, cg_crd, avg_crd, cb_crd] 
  #                 sc_vector]  

  # Create DFG object for MPI
  mpi = multiprocessing.Pool(processes = multiprocessing.cpu_count())
  Ref = DFGTorsions(Ref_Coords)
  Tmp = mpi.map(DFGTorsions, Tgt_Coords)
  mpi.close()
  mpi.join()

  Tgt_List = [x for x in Tmp if x is not None]
  print('\n ## DFG-motif Vector return: {0}\n'.format(len(Tgt_List)))
  ExportDFGMeasure(Ref, Tgt_List, Data, output)


##########################################################################
## Compare template and model DFG-motif cross-product vectors
def ExportDFGMeasure( Ref, Tgt_List, Data, output ):

#  fo = open(output+'.dfg_vec.txt', 'w')
#  fo.write('#<DFG> PDB_h | Resi  | p1.p1x | p2.p2x | DFG-type\n')

  # Reference DFG parameters
  pdb_id_x, resi_id_x, p1x, p2x, d0_presx, d1_presx = Ref
  state = DFGState( np.vdot(p1x,p1x), np.vdot(p2x,p2x) )
  Ref_Ref = [ pdb_id_x, resi_id_x, np.vdot(p1x,p1x), np.vdot(p2x,p2x),
              state, d0_presx, d1_presx ]

  Tgt_Tmp = []
  for Tgt in Tgt_List:
    pdb_id, resi_id, p1, p2, d0_pres, d1_pres = Tgt
    if p1 is None or p2 is None:
      p1p1x, p2p2x = None, None
    else:
      p1p1x = np.vdot(p1, p1x)
      p2p2x = np.vdot(p2, p2x)
    state = DFGState(p1p1x, p2p2x)
    Tgt_Tmp.append([ pdb_id, resi_id, p1p1x, p2p2x, state, d0_pres, d1_pres ])
  Tgt_Inp = sorted(Tgt_Tmp, key=lambda x: (x[4], x[0]))

  # Add Ref to the target list for printing
  All_Data = [Ref_Ref]
  All_Data.extend(Tgt_Inp)

  for idx, V in enumerate(All_Data):
    pdb_id, resi_id, p1p1x, p2p2x, state, d0_x, d1_x =  V

    Data[pdb_id]['p1p1x']  = p1p1x
    Data[pdb_id]['p2p2x']  = p2p2x
    Data[pdb_id]['dfg_st'] = state
    Data[pdb_id]['d0_x']   = d0_x
    Data[pdb_id]['d1_x']   = d1_x

#    dihe  = np.arccos( np.vdot(p1/VecMag(p1), p2/VecMag(p2)) )
#    print('Mod Dihe {0} -- {1}'.format(pdb_id, dihe))
#    if idx > 0: key = 'VEC'
#    else: key = 'REF'
#    if   p1p1x is None and p2p2x is not None:
#      fo.write('<{0}> {1:6} | {2:5} |  None  | {3:6.3f} | {4}\n'.format(
#                key, pdb_id, resi_id, p2p2x, state ))
#    elif p1p1x is not None and p2p2x is None:
#      fo.write('<{0}> {1:6} | {2:5} | {3:6.3f} |  None  | {4}\n'.format(
#                key, pdb_id, resi_id, p1p1x, state ))
#    elif p1p1x is None and p2p2x is None:
#      fo.write('<{0}> {1:6} | {2:5} |  None  |  None  | {3}\n'.format(
#                key, pdb_id, resi_id, state ))
#    else:
#      fo.write('<{0}> {1:6} | {2:5} | {3:6.3f} | {4:6.3f} | {5}\n'.format( 
#                key, pdb_id, resi_id, p1p1x, p2p2x, state ))
#  fo.close()


##########################################################################
##
def DFGState( p1p1x, p2p2x ):

  ## Model PDB has same DFG- config as template DFG-in:     'in'
  ## Model PDB has opposite DFG- config as template DFG-in: 'out'
  ## Model PDB has undefined DFG- config:                   'random'
  if p1p1x is None or p2p2x is None:
    return 'missing DFG'
  elif p1p1x > -0.005 and p2p2x > 0.00:
    return 'in'
  elif p1p1x < -0.125 and p2p2x < -0.05:
    return 'out-1'
  elif p1p1x > -0.125 and p1p1x < -0.0005 and p2p2x < -0.05:
    return 'out-2'
  elif p1p1x > -0.0005 and p1p1x < 0.005 and p2p2x < 0.00:
    return 'F-flip'
  else:
    return 'random'


##########################################################################
##
def DFGTorsions( Input ):
  # Input_Coords = [pdb_name, H_Crds, N_Crds, C_Crds, G_Crds, R_Crds, T_Crds]
  #     x_Coords = [resname, resid, bb_crds, ca_crd, cg_crd, avg_crd, cb_crd] 
  pdb_name, Pre_Coords = Input[0], Input[3]

  # Check for missing residues
  if Pre_Coords is None:
    print('  # DFG Warning: No DFG residue available: {0}'.format(pdb_name))
    return None
  for idx, Seq in enumerate(Pre_Coords):
    if Seq is None:
      print('  # DFG Warning: Missing resid: '+pdb_name+' '+str(idx+1))
      return None

  # reformat the data array
  CA_Coords = np.asarray(zip(*Pre_Coords)[3])
  CG_Coords = np.asarray(zip(*Pre_Coords)[4])
  CB_Coords = np.asarray(zip(*Pre_Coords)[6])

  center = ArrayCent(len(CA_Coords))
  res_id = AA(Pre_Coords[center][0])+str(Pre_Coords[center][1])

  # Check number of residue collected for DFG-motif. If D+1 residue is not
  # collected (only D's N,CA,C atoms, no D+1's N,CA,C atoms), report and skip
  d0_x, d1_x = 1, 1
  if len(Pre_Coords) < 3:
    print('  # DFG Warning: Missing D+1 residue for measurement: '+pdb_name)
    p1, p2, d0_x, d1_x = None, None, None, None
  else:
    ## Take Asp 'D' and Phe 'D+1' from DFG to measure; 
    ## if CG not available, use CB as substitute but mark as no CG presence
    ## r1 = Asp 'CG', r2 = Asp 'CA'
    ## r3 = Phe 'CA', r4 = Phe 'CG'
    if CG_Coords[center] is None:
      p0_x = None
      print('  # DFG Warning: Substitute [D] CG with CB: '+pdb_name)
      CG_Coords[center] = CB_Coords[center]

    if CG_Coords[center+1] is None:
      d1_x = None
      print('  # DFG Warning: Substitute [D+1] CG with CB: '+pdb_name)
      CG_Coords[center+1] = CB_Coords[center+1]

    p1, p2 = CalculateVector( CG_Coords[center],   CA_Coords[center],
                              CA_Coords[center+1], CG_Coords[center+1],
                              pdb_name )
  
  return [ pdb_name, res_id, p1, p2, d0_x, d1_x ]


##########################################################################
## Take in coordinates, calculate vectors among the coords, generate
## Cross-Products of the pairs
## Asp-CG (r1), Asp-CA (r2), Phe-CA (r3), Phe-CG (r4)
def CalculateVector(r1, r2, r3, r4, pdb_name):
  
  try:
    r21 = r1-r2   # (AspCG-AspCA)
  except TypeError:
    print('  # DFG Warning: No [D] CG/CB for AspCG-AspCA: '+pdb_name)
    return None, None
  r23 = r3-r2   # (AspCA-PheCA)
  r32 = r2-r3   # (PheCA-AspCA)
  try:
    r34 = r4-r3   # (PheCG-PheCA)
  except TypeError:
    print('  # DFG Warning: No [D+1] CG/CB for PheCG-PheCA: '+pdb_name)
    return None, None

  p1 = np.cross(r21,r23)/(VecMag(r21)*VecMag(r23))
  p2 = np.cross(r34,r32)/(VecMag(r34)*VecMag(r32))
  return p1, p2


###########################################################################
## Extract the Phi/Psi torsional angles of the residues from BioPython
def Phi_Psi(chain_obj):
  PP_List = []
  for chain in chain_obj:
    polypep = PPBuilder().build_peptides(chain)
    for poly in polypep:
      phi_psi = poly.get_phi_psi_list()
      for pp in phi_psi:
        Set = []
        for num in list(pp):
          if num is not None:
            Set.append(np.rint(np.degrees(num)).astype(int))
          else:
            Set.append(0)
        PP_List.append(Set)
  return PP_List


###########################################################################
##
##  Peter M.U. Ung @ MSSM
##
##  v0.1    17.01.30
##
##  v1.0    17.02.02
##  v1.1    17.02.05
##  v2.0    17.02.21	add data collection object
##  v3.0    17.05.08    use CB for DFG vector if CG coords not available
##
