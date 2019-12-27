#!/usr/bin/python


import numpy as np
from pathos import multiprocessing
from tqdm import tqdm
from aa_residue    import *
from x_helix_axis  import *
from Bio.PDB.Polypeptide import PPBuilder
np.seterr(invalid='ignore')

##########################################################################
##
def DFGTorsionAngle( Ref_Coords, Tgt_Coords, Data, output ):

  # Input_Coords = [pdb_name, H_Crds, N_Crds, C_Crds, G_Crds, R_Crds, T_Crds]
  #     x_Coords = [resname, resid, bb_crds, ca_crd, cg_crd, avg_crd, cb_crd] 
  #                 sc_vector]  

  print('##################################################################\n')

  # Create DFG object for MPI
  mpi = multiprocessing.Pool(processes = multiprocessing.cpu_count())
  Ref = DFGTorsions(Ref_Coords)
  print(Ref)
#  Tmp = [DFGTorsions(Tgt) for Tgt in Tgt_Coords]
  Tmp = [x for x in tqdm(mpi.imap_unordered(DFGTorsions, Tgt_Coords),total=len(Tgt_Coords))]
  mpi.close()
  mpi.join()

  Tgt_List = [x for x in Tmp if x is not None]
  print('\n ## DFG-motif Vector return: {0}\n'.format(len(Tgt_List)))
  ExportDFGMeasure(Ref, Tgt_List, Data, output)


##########################################################################
## Compare template and model DFG-motif cross-product vectors
def ExportDFGMeasure( Ref, Tgt_List, Data, output ):

  # Reference DFG parameters
  pdb_id_x, resi_id_x, p1x, p2x, v1x, v2x, r3x, d0_presx, d1_presx = Ref
  p_state = DFGState( np.vdot(p1x,p1x), np.vdot(p2x,p2x) )
  v_state = DFGState( np.vdot(v1x,v1x), np.vdot(v2x,v2x) )

  Ref_Ref = [ pdb_id_x, resi_id_x, np.vdot(p1x,p1x), np.vdot(p2x,p2x),
              p_state, np.vdot(v1x,v1x), np.vdot(v2x,v2x), v_state,
              r3x, d0_presx, d1_presx ]

  Tgt_Tmp = []
  for Tgt in Tgt_List:
    pdb_id, resi_id, p1, p2, v1, v2, r3, d0_pres, d1_pres = Tgt
    if p1 is None or p2 is None:	# for cross product DFG vectors
      p1p1x, p2p2x = None, None
    else:
      p1p1x = np.vdot(p1, p1x)
      p2p2x = np.vdot(p2, p2x)
    p_state = DFGState(p1p1x, p2p2x)

    if v1 is None or v2 is None:	# for DFG vectors
      v1v1x, v2v2x = None, None
    else:
      v1v1x = np.vdot(v1, v1x)
      v2v2x = np.vdot(v2, v2x)
    v_state = DFGState(v1v1x, v2v2x)

    if r3 is None:
      r3r3x = None
    else:
      r3r3x = np.vdot(r3, r3x)
  
    Tgt_Tmp.append([ pdb_id, resi_id, p1p1x, p2p2x, p_state,
                     v1v1x, v2v2x, v_state, r3r3x, d0_pres, d1_pres ])
  Tgt_Inp = sorted(Tgt_Tmp, key=lambda x: (x[4], x[0]))	# sort on crossed vec

  # Add Ref to the target list for printing
  All_Data = [Ref_Ref]
  All_Data.extend(Tgt_Inp)

  for idx, V in enumerate(All_Data):
    pdb_id, resi_id, p1p1x, p2p2x, p_state, v1v1x, v2v2x, v_state, r3r3x, d0_x, d1_x =  V

    Data[pdb_id]['p1p1x']  = p1p1x
    Data[pdb_id]['p2p2x']  = p2p2x
    Data[pdb_id]['v1v1x']  = v1v1x
    Data[pdb_id]['v2v2x']  = v2v2x
    Data[pdb_id]['r3r3x']  = r3r3x
    Data[pdb_id]['dfg_vs'] = v_state
    Data[pdb_id]['dfg_st'] = p_state
    Data[pdb_id]['d0_x']   = d0_x
    Data[pdb_id]['d1_x']   = d1_x


##########################################################################
##
def DFGTorsions( Input ):
  # Input_Coords = [pdb_name, H_Crds, N_Crds, C_Crds, G_Crds, F_Crds, T_Crds]
  #     x_Coords = [resname, resid, bb_crds, ca_crd, cg_crd, avg_crd, cb_crd] 
  pdb_name, D_Coords, F_Coords = Input[0], Input[3], Input[5]

  # Check for missing DFG-D or DFG-F residue
  if D_Coords is None:
    print('\n  #2# DFG Warning: No DFG-D residue available: {0}'.format(pdb_name))
    return None
  if F_Coords is None:
    print('\n  #2# DFG Warning: No DFG-F residue available: {0}'.format(pdb_name))
    return None

  for idx, Seq in enumerate(D_Coords):
    if Seq is None:
      print('\n  #2# DFG Warning: Missing DFG-D resid: '+pdb_name+' '+str(idx+1))
      return None
  for idx, Seq in enumerate(F_Coords):
    if Seq is None:
      print('\n  #2# DFG Warning: Missing DFG-F resid: '+pdb_name+' '+str(idx+1))
      return None


  # reformat the data array for DFG-D and DFG-F
  D_CA_Coords = np.asarray(zip(*D_Coords)[3])
  D_CG_Coords = np.asarray(zip(*D_Coords)[4])
  D_CB_Coords = np.asarray(zip(*D_Coords)[6])
  center = ArrayCent(len(D_CA_Coords))
  res_id = AA(D_Coords[center][0])+str(D_Coords[center][1])

  F_CA_Coords = np.asarray(zip(*F_Coords)[3])
  F_CG_Coords = np.asarray(zip(*F_Coords)[4])
  F_CB_Coords = np.asarray(zip(*F_Coords)[6])
  f_center = ArrayCent(len(F_CA_Coords))
  f_res_id = AA(F_Coords[f_center][0])+str(F_Coords[f_center][1])


  ## Take Asp 'D' and Phe 'D+1' from DFG to measure; 
  ## if CG not available, use CB as substitute but mark as no CG presence
  ## r1 = Asp 'CG', r2 = Asp 'CA'
  ## r3 = Phe 'CA', r4 = Phe 'CG'
  d0_x, d1_x = 1, 1

  if D_CG_Coords[center] is None:
    p0_x = None
    print('\n  #1# DFG Warning: Substitute [D] CG with CB: '+pdb_name)
    D_CG_Coords[center] = D_CB_Coords[center]

  if F_CG_Coords[f_center] is None:
    d1_x = None
    print('\n  #1# DFG Warning: Substitute [D+1] CG with CB: '+pdb_name)
    F_CG_Coords[f_center] = F_CB_Coords[f_center]

  p1, p2, v1, v2, r3  = CalculateVector( 
                            D_CG_Coords[center],   D_CA_Coords[center],
                            F_CA_Coords[f_center], F_CG_Coords[f_center],
                            pdb_name )
  
  return [ pdb_name, res_id, p1, p2, v1, v2, r3, d0_x, d1_x ]


##########################################################################
## Take in coordinates, calculate vectors among the coords, generate
## Cross-Products of the pairs
## Asp-CG (r1), Asp-CA (r2), Phe-CA (r3), Phe-CG (r4)
def CalculateVector(r1, r2, r3, r4, pdb_name):
  
  try:
    r21 = np.array(r1-r2, dtype=np.float64)   # (AspCG-AspCA)
  except TypeError:
    print('\n  #2# DFG Warning: No [D] CG/CB for AspCG-AspCA: '+pdb_name)
    return [None, None, None, None, None]
  r23 = np.array(r3-r2, dtype=np.float64)   # (AspCA-PheCA)
  r32 = np.array(r2-r3, dtype=np.float64)   # (PheCA-AspCA)
  try:
    r34 = np.array(r4-r3, dtype=np.float64)   # (PheCG-PheCA)
  except TypeError:
    print('\n  #2# DFG Warning: No [D+1] CG/CB for PheCG-PheCA: '+pdb_name)
    return [None, None, None, None, None]

  # V = vector of the side chain, P = cross product vector of sidechain/D:D+1
  p1 = np.cross(r21/VecMag(r21),r23/VecMag(r23))
  p2 = np.cross(r34/VecMag(r34),r32/VecMag(r32))
  v1 = r21/VecMag(r21)
  v2 = r34/VecMag(r34)
  r3 = r23/VecMag(r23)

  return [ p1/VecMag(p1), p2/VecMag(p2), v1, v2, r3 ]


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
def DFGState( p1p1x, p2p2x ):

  ## Model PDB has same DFG- config as template DFG-in:     'in'
  ## Model PDB has opposite DFG- config as template DFG-in: 'out'
  ## Model PDB has undefined DFG- config:                   'random'
  if p1p1x is None or p2p2x is None:
    return 'missing DFG'
  elif p1p1x > 0.005 and p2p2x > 0.05:
    return 'in'
  elif p1p1x < -0.125 and p2p2x < -0.125:
    return 'out'
  else:
    return 'random'


##########################################################################
##
##  Peter M.U. Ung @ MSSM
##
##  v0.1    17.01.30
##
##  v1.0    17.02.02
##  v1.1    17.02.05
##  v2.0    17.02.21	add data collection object
##  v3.0    17.05.08    use CB for DFG vector if CG coords not available
##  v4.0    17.10.07	add non-crossed vectors of D/D+1
##  v5.0    18.03.15  update to include DFG-F detection and correction
##
