#!/usr/bin/env python3

from pathos import multiprocessing
from tqdm import tqdm
from x_dfg_torsion import *
from x_helix_axis  import *
from x_fasta_parse import *
from x_pdb_extract import *

##########################################################################
## Use MPI to run calculation of spherical angle and distance calculations 
## for H-helix, N-domain, and C-domain. Results are exported into files
def DomainDistances( Ref_Coords, PDB_Coords, RefReg2, Reg2, Data, output ):

  # Input_Coords = [pdb_name, H_Crds, N_Crds, C_Crds, G_Crds, R_Crds, T_Crds]
  #     x_Coords = [resname, resid, bb_crds, ca_crd, cg_crd, avg_crd, cb_crd] 

  print('##################################################################\n')

  # Create distance object for MPI
  mpi = multiprocessing.Pool(processes = multiprocessing.cpu_count())
  Ref = CalculateDist( [Ref_Coords, RefReg2] )
  Tmp = [CalculateDist([Tgt,Reg2[idx]]) for idx,Tgt in enumerate(PDB_Coords)]
#  Tmp = [x for x in tqdm(mpi.imap_unordered( CalculateDist, list(zip(PDB_Coords, Reg2))),
#                                       total=len(Reg2))]
  mpi.close()
  mpi.join()

  Tgt_List = [x for x in Tmp if x is not None]
  print('\n ## Domain Distances return: {0}\n'.format(len(Tgt_List)))
  CollectDomain(Ref, Tgt_List, Data)


##########################################################################
## Calculate C-helix-domains parameters. For each element, use the center
## residue in the array. Since input is odd number and 3 atoms for each 
## residue, array center be should 'CA' of center residue
## e.g.>  HTLN|E|KRIL
def CalculateDist( Input ):

  None_Coord = [None, None, None, None, None, None, None, None, None,
                None, None, None, None, None, None ]

  Coords, Reg2 = Input
  pdb_id, Pre_Coords = Coords[0:2]

  if Pre_Coords is None:
    print('\n  #2# Domain Warning: Too short to calculate: {0}'.format(pdb_id))
    return None_Coord
  for idx, Seq in enumerate(Pre_Coords):
    if Seq is None:
      print('\n  #2# Domain Warning: Missing residue: {0}\t{1}'.format(pdb_id,idx+1))
      return None_Coord

  if RetreiveCoords(Coords) is not None:
    H_ca, H_cb, H_cg, H_cd, N_ca, N_cg, N_cd, C_ca, C_cg = RetreiveCoords(Coords)
  else:
    print('\n  #2# Domain Warning: Cannot retreive coordinates: '+pdb_id)
    return None_Coord

  ## Calculate CA-CA distance between N/C-domains (Lys/Asp) to C-helix Glu
  # Calculate CG-CG distance fo  dist_N_H may be better - 17.10.09
  # Calculate CD-CD distnace for dist_N_H should be ok  - 17.10.13
  dist_C_H = Distance(H_ca, C_ca)
  dist_N_C = Distance(N_ca, C_ca)
  dist_N_H = Distance(H_cg, N_cg)
#  dist_N_H = Distance(H_cd, N_cd)
#  dist_N_H = Distance(H_ca, N_ca)

################
  ## Calculate spherical angles from N/C-domain to origin
  n_phi, n_psi = SphericalAngles(N_ca)
  c_phi, c_psi = SphericalAngles(C_ca)

################
  ## Calculate N(C)-lobe/C-helix center/C-Glu sc angle 

  # Helix center corresponding to center residue
  try:
    H_rg = Reg2[ArrayCent(len(Reg2))]
  except TypeError:
    print('\n  #2#  Domain Warning: Reg2 (C-helix linear regression) failed: '+pdb_id)
    H_rg = None

  # Check if side chain is there by using H_cb presence
  if H_cb is None or H_rg is None:
    print('\n  #2# Domain Warning: Helix data is missing: '+Coords[0]+' '+Coords[1][ArrayCent(len(Coords[1]))][0]+str(Coords[1][ArrayCent(len(Coords[1]))][1]))
    angl_NHs, angl_CHs = None, None
  else:
    if H_cd is None:
      print('\n  #1# Domain -- H_CD not available, use H_CG for angle: '+pdb_id)
      H_cd = H_cg
    elif H_cg is None:
      print('\n  #1# Domain -- H_CD/H_CG not available, use H_CB for angle: '+pdb_id)
      H_cd = H_cb
    # Calculate angle between N-Lys_CA/C-Asp_CA -- H-Glu_CA -- H-Glu_CD
    angl_NHs = VectorAngle( (H_rg - N_ca), (H_rg - H_cg) )
    angl_CHs = VectorAngle( (H_rg - C_ca), (H_rg - H_cg) )

###############
  # Gatekeeper residue
  G_Crds = Coords[4][ArrayCent(len(Coords[4]))]
  g_resi = '{0} {1}'.format(G_Crds[0], G_Crds[1])
  g_coord= G_Crds[5]

  return [ Coords[0], H_ca, H_cg, N_ca, C_ca,
                      n_phi, n_psi, c_phi, c_psi,
                      dist_N_H, dist_N_C, dist_C_H,
                      angl_NHs, angl_CHs,
                      g_resi, g_coord ]


##########################################################################
def CollectDomain(Ref, Tgt_List, Data):
  Tgt_List.append(Ref)
  for Tgt in Tgt_List:
    if Tgt is None:
      continue
    if Tgt[0] is None:
      continue

    if Tgt[2] is not None:
      Data[Tgt[0]]['h_cg']    = Distance(Tgt[2],Ref[2])

#    Data[Tgt[0]]['h_ca']    = Distance(Tgt[1],Ref[1])
#    Data[Tgt[0]]['n_ca']    = Distance(Tgt[3],Ref[3])
#    Data[Tgt[0]]['c_ca']    = Distance(Tgt[4],Ref[4])
#    Data[Tgt[0]]['n_phi']   = Tgt[5]
#    Data[Tgt[0]]['n_psi']   = Tgt[6]
#    Data[Tgt[0]]['c_phi']   = Tgt[7]
#    Data[Tgt[0]]['c_psi']   = Tgt[8]
    Data[Tgt[0]]['dist_NH'] = Tgt[9]
    Data[Tgt[0]]['dist_NC'] = Tgt[10]
    Data[Tgt[0]]['dist_CH'] = Tgt[11]
    Data[Tgt[0]]['ang_NHs'] = Tgt[12]
    Data[Tgt[0]]['ang_CHs'] = Tgt[13]
    Data[Tgt[0]]['g_resi']  = Tgt[14]
#    Data[Tgt[0]]['g_coord'] = Tgt[15]


##########################################################################
def RetreiveCoords( Coords ):
  # Input_Coords = [pdb_name, H_Crds(1), N_Crds(2), C_Crds(3), G_Crds(4),
  #                           F_Crds(5), T_Crds(6)]
  #     x_Coords = [resname, resid, bb_crds(2), ca_crd(3), cg_crd(4), 
  #                                 avg_crd(5), cb_crd(6), cd_crd(7) ]  

  pdb_id = Coords[0]

  # Check for missing residues
  for idx, Seq in enumerate(Coords):
    name = ['0', 'Helix', 'N-lobe', 'C-lobe','Gate', 'DFG-F', 'C-spine']
    if Seq is None:
      print('\n  #2# Domain Warning: Missing '+name[idx]+' residue: '+pdb_id)
      return None

  # Extract Helix, C-dom, N-dom CA coordinates of the center residue
  # if CA atom is not available, use CA atom as substitute
  if Coords[1][ArrayCent(len(Coords[1]))][3] is not None:
    H_ca = Coords[1][ArrayCent(len(Coords[1]))][3]
  else:
    print('\n  #1# Domain -- C-Helix residue missing CA: '+pdb_id)
  if Coords[1][ArrayCent(len(Coords[1]))][6] is not None:
    H_cb = Coords[1][ArrayCent(len(Coords[1]))][6]
  else:
    print('\n  #1# Domain -- C-Helix residue missing CA: '+pdb_id)
    H_cb = Coords[1][ArrayCent(len(Coords[1]))][3]

  if Coords[1][ArrayCent(len(Coords[1]))][4] is not None:
    H_cg = Coords[1][ArrayCent(len(Coords[1]))][4]
  elif Coords[1][ArrayCent(len(Coords[1]))][6] is not None:
    print('\n  #1# Domain -- C-Helix residue missing CG, use CB: '+pdb_id)
    H_cg = Coords[1][ArrayCent(len(Coords[1]))][6]
  else:
    print('\n  #1# Domain -- C-Helix residue missing CG, CB, use CA: '+pdb_id)
    H_cg = Coords[1][ArrayCent(len(Coords[1]))][3]

  if Coords[1][ArrayCent(len(Coords[1]))][7] is not None:
    H_cd = Coords[1][ArrayCent(len(Coords[1]))][7]
  elif Coords[1][ArrayCent(len(Coords[1]))][4] is not None:
    print('\n  #1# Domain -- C-Helix residue missing CD, use CG: '+pdb_id)
    H_cd = Coords[1][ArrayCent(len(Coords[1]))][4]
  elif Coords[1][ArrayCent(len(Coords[1]))][6] is not None:
    print('\n  #1# Domain -- C-Helix residue missing CD,CG, use CB: '+pdb_id)
    H_cd = Coords[1][ArrayCent(len(Coords[1]))][6]
  else:
    print('\n  #1# Domain -- C-Helix residue missing CD,CG,CB use CA: '+pdb_id)
    H_cd = Coords[1][ArrayCent(len(Coords[1]))][3]

  if Coords[2][ArrayCent(len(Coords[2]))][3] is not None:
    N_ca = Coords[2][ArrayCent(len(Coords[2]))][3]
  else:
    print('\n  #1# Domain -- N-lobe residue missing CA: '+pdb_id)
  if Coords[2][ArrayCent(len(Coords[2]))][4] is not None:
    N_cg = Coords[2][ArrayCent(len(Coords[2]))][4]
  else:
    print('\n  #1# Domain -- N-lobe residue missing CG: '+pdb_id)
    N_cg = Coords[2][ArrayCent(len(Coords[2]))][3]

  if Coords[1][ArrayCent(len(Coords[1]))][7] is not None:
    N_cd = Coords[1][ArrayCent(len(Coords[1]))][7]
  elif Coords[1][ArrayCent(len(Coords[1]))][4] is not None:
    print('\n  #1# Domain -- N-lobe residue missing CD, use CG: '+pdb_id)
    N_cd = Coords[1][ArrayCent(len(Coords[1]))][4]
  elif Coords[1][ArrayCent(len(Coords[1]))][6] is not None:
    print('\n  #1# Domain -- N-lobe residue missing CD,CG, use CB: '+pdb_id)
    N_cd = Coords[1][ArrayCent(len(Coords[1]))][6]
  else:
    print('\n  #1# Domain -- N-lobe residue missing CD,CG,CB use CA: '+pdb_id)
    N_cd = Coords[1][ArrayCent(len(Coords[1]))][3]

  if Coords[3][ArrayCent(len(Coords[3]))][3] is not None:
    C_ca = Coords[3][ArrayCent(len(Coords[3]))][3]
  else:
    print('\n  #1# Domain -- C-lobe residue missing CA: '+pdb_id)
  if Coords[3][ArrayCent(len(Coords[3]))][4] is not None:
    C_cg = Coords[3][ArrayCent(len(Coords[3]))][4]
  else:
    print('\n  #1# Domain -- C-lobe residue missing CG: '+pdb_id)
    C_cg = Coords[3][ArrayCent(len(Coords[3]))][3]

  return [ H_ca, H_cb, H_cg, H_cd, N_ca, N_cg, N_cd, C_ca, C_cg ]


##########################################################################
#
#   Peter M.U. Ung  @ MSSM
#   
#   v0.1    17.01.30
#   v0.2    17.02.01
#   
#   v1.0    17.02.02
#   v2.0    17.02.21	- add data collection object
#   v3.0    17.03.10
#   v4.0    17.10.09	dist_NH change from CA-CA to sc-sc
#   v5.0    17.12.21    only calculate the relevant parameters
#   v6.0    18.03.12    fixed bugs
#
