#!/usr/bin/python

import re,os,glob,sys
import numpy as np
import multiprocessing
from x_dfg_torsion import *
from x_helix_axis  import *
from x_fasta_parse import *
from x_pdb_extract import *
from CommonUtility import *


##########################################################################
## Use MPI to run calculation of spherical angle and distance calculations 
## for H-helix, N-domain, and C-domain. Results are exported into files
def DomainDistances( Ref_Coords, PDB_Coords, RefReg2, Reg2, Data, output ):

  # Input_Coords = [pdb_name, H_Crds, N_Crds, C_Crds, G_Crds, R_Crds, T_Crds]
  #     x_Coords = [resname, resid, bb_crds, ca_crd, cg_crd, avg_crd, cb_crd] 

  # Create distance object for MPI
  mpi = multiprocessing.Pool(processes = multiprocessing.cpu_count())
  Ref = CalculateDist( [Ref_Coords, RefReg2] )
#  Tmp = [CalculateDist( [Tgt, Reg2[idx]] ) for idx, Tgt in enumerate(PDB_Coords)]
  Tmp = mpi.map( CalculateDist, zip(PDB_Coords, Reg2) )
  mpi.close()
  mpi.join()

  Tgt_List = [x for x in Tmp if x is not None]
  print('\n ## Domain Distances return: {0}\n'.format(len(Tgt_List)))
#  ExportDistMeasure(Ref, Tgt_List, output)
  CollectDomain(Ref, Tgt_List, Data)


##########################################################################
## Calculate C-helix-domains parameters. For each element, use the center
## residue in the array. Since input is odd number and 3 atoms for each 
## residue, array center be should 'CA' of center residue
## e.g.>  HTLN|E|KRIL
def CalculateDist( Input ):
  Coords, Reg2 = Input
  pdb_id = Coords[0]

  # Input_Coords = [pdb_name, H_Crds(1), N_Crds(2), C_Crds(3), G_Crds(4),
  #                           R_Crds(5), T_Crds(6)]
  #     x_Coords = [resname, resid, bb_crds(2), ca_crd(3), cg_crd(4), 
  #                                 avg_crd(5), cb_crd(6)]  

  # Check for missing residues
  for Seq in Coords:
    if Seq is None:
      print('  # Domain Warning: Missing residue: '+Coords[0])
      return None

  # Extract Helix, C-dom, N-dom CA coordinates of the center residue
  H_ca = Coords[1][ArrayCent(len(Coords[1]))][3]
  C_ca = Coords[3][ArrayCent(len(Coords[3]))][3]
  try:
    N_ca = Coords[2][ArrayCent(len(Coords[2]))][3]
  except TypeError:
    print(pdb_id)
    

  # Calculate spherical angles from N/C-domain to origin
  # Calculate CA-CA distance between N/C-domains (Lys/Asp) to C-helix Glu
  n_phi, n_psi = SphericalAngles(N_ca)
  c_phi, c_psi = SphericalAngles(C_ca)
  dist_N_H = Distance(H_ca, N_ca)
  dist_C_H = Distance(H_ca, C_ca)
  dist_N_C = Distance(N_ca, C_ca)

################
  # Helix center corresponding to center residue
  H_rg = Reg2[ArrayCent(len(Reg2))]

  # Extract the CG coordinates of the side chain of residue; if CG is missing,
  # use CB to substitute; CB is in all AA except Gly
  H_cg = Coords[1][ArrayCent(len(Coords[1]))][4]
  H_cb = Coords[1][ArrayCent(len(Coords[1]))][6]
  if H_cb is None:
    print('  # Domain Warning: Helix side chain missing: '+Coords[0]+' '+Coords[1][ArrayCent(len(Coords[1]))][0]+str(Coords[1][ArrayCent(len(Coords[1]))][1]))
    angl_NHs, angl_CHs = None, None
  else:
    if H_cg is None:
      print('  # Domain Warning: H_CG not available: '+pdb_id)
      H_cg = H_cb
    # Calculate angle between N-Lys_CA/C-Asp_CA -- H-Glu_CA -- H-Glu_CG
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

    if Tgt[2] is not None:
      Data[Tgt[0]]['h_cg']    = Distance(Tgt[2],Ref[2])

    Data[Tgt[0]]['h_ca']    = Distance(Tgt[1],Ref[1])
    Data[Tgt[0]]['n_ca']    = Distance(Tgt[3],Ref[3])
    Data[Tgt[0]]['c_ca']    = Distance(Tgt[4],Ref[4])
    Data[Tgt[0]]['n_phi']   = Tgt[5]
    Data[Tgt[0]]['n_psi']   = Tgt[6]
    Data[Tgt[0]]['c_phi']   = Tgt[7]
    Data[Tgt[0]]['c_psi']   = Tgt[8]
    Data[Tgt[0]]['dist_NH'] = Tgt[9]
    Data[Tgt[0]]['dist_NC'] = Tgt[10]
    Data[Tgt[0]]['dist_CH'] = Tgt[11]
    Data[Tgt[0]]['ang_NHs'] = Tgt[12]
    Data[Tgt[0]]['ang_CHs'] = Tgt[13]
    Data[Tgt[0]]['g_resi']  = Tgt[14]
    Data[Tgt[0]]['g_coord'] = Tgt[15]


##########################################################################
## Export spherical angle and distance results into files
## All numbers are relative to a reference structure (1atp)
def ExportDistMeasure( Ref, Tgt_List, output ):
#          [ pdb_name, h_ca(1),      h_cg(2),     n_ca(3),     c_ca(4),
#                      n_phi(5),     n_psi(6),    c_phi(7),    c_psi(8),
#                      dist_N_H(9),  dist_N_C(10),dist_C_H(11),
#                      angl_NHs(12), angl_CHs(13),
#                      g_resi(14),   g_coord(15) ]

  # Write out relative atomic displacement
  cg   = open(output+'.cg_dist.txt', 'w')
  cg.write('## Relative side chain distance to reference\n')
  cg.write('## Reference {0}: Helix CA | Helix CG | N-Dom CA | C-Dom CA\n'.format( Ref[0] ))

  # Write out relative spherical angles for domains
  angl = open(output+'.dom-angle.txt', 'w')
  angl.write('## Spherical Angle of N-dom, C-dom\n')
  angl.write('## Reference {0} N(Phi:Psi): {1:5.1f}  {2:5.1f}  C(Phi:Psi): {3:5.1f}  {4:5.1f}\n'.format(
             Ref[0], Ref[5], Ref[6], Ref[7], Ref[8] ))

  # Write out relative inter-lobe distances
  dist = open(output+'.dom-dist.txt', 'w')
  dist.write('## Dist between N-dom_C-heix, C-dom_C-helix, N-dom_C-dom\n')
  dist.write('## Reference {0} N-H: {1:5.1f}  C-H: {2:5.1f} N-C: {3:5.1f}\n'.format(        Ref[0], Ref[9], Ref[11], Ref[10] ))

  # Write out relative C-helix-Ndom/Cdom angles
  helx = open(output+'.h_dom-angle.txt', 'w')
  helx.write('## Angle between N-Lys_H-Cent_H-Glu, C-Asp_H-Cent_H-Glu\n')
  helx.write('## Reference {0} N-Lys_H-Cent_H-Glu: {1:5.1f}  C-Asp_H-Cent_H-Glu: {2:5.1f}\n'.format(Ref[0], Ref[12], Ref[13] ))


  for Tgt in Tgt_List:
    cg.write('{0}: {1:5.1f}  {2:5.1f}  {3:5.1f}  {4:5.1f}\n'.format(
        Tgt[0], VecMag(Tgt[1]-Ref[1]), VecMag(Tgt[2]-Ref[2]), 
        VecMag(Tgt[3]-Ref[3]), VecMag(Tgt[4]-Ref[4]) ) )

    angl.write('{0} N(Phi:Psi): {1:5.1f}  {2:5.1f}  C(Phi:Psi): {3:5.1f}  {4:5.1f}\t|\t{5:5.1f}  {6:5.1f}  {7:5.1f}  {8:5.1f}\n'.format( 
        Tgt[0], Tgt[5], Tgt[6], Tgt[7], Tgt[8],
          (Tgt[5]-Ref[5]), (Tgt[6]-Ref[6]), (Tgt[7]-Ref[7]), (Tgt[8]-Ref[8]) ))

    dist.write('{0} (N-H:C-H:N-C): {1:5.1f}  {2:5.1f}  {3:5.1f}\t|\t{4:5.1f}  {5:5.1f}  {6:5.1f}\n'.format(
        Tgt[0], Tgt[9], Tgt[11], Tgt[10],
          (Tgt[9]-Ref[9]), (Tgt[11]-Ref[11]), (Tgt[10]-Ref[10]) ))

    if Tgt[12] is None:
      helx.write('{0} (N-H:H-s|C-H:H-s): {1:5.1f}  \t|\t{3:5.1f}  {4:5.1f}\n'.format(
        Tgt[0]))
    else:
      helx.write('{0} (N-H:H-s|C-H:H-s): {1:5.1f}  {2:5.1f}\t|\t{3:5.1f}  {4:5.1f}\n'.format(
        Tgt[0], Tgt[12], Tgt[13], 
          (Tgt[12]-Ref[12]), (Tgt[13]-Ref[13]) ))

  angl.close()
  dist.close()
  helx.close()
  cg.close()

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
