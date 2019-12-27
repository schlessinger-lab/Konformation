#!/usr/bin/python

import re,os,glob,sys
import numpy as np
import multiprocessing

from x_helix_axis  import *


##########################################################################
class RCSpines( object ):
  
  def __init__( self, ref_r=None, ref_c=None ):
    self.ref_r = ref_r
    self.ref_c = ref_c

  def __call__(self, Input):
    return self.spine_measurement(Input)

#################
  def spine_measurement( self, Input ):

  # Input_Coords = [pdb_name, H_Crds, N_Crds, C_Crds, G_Crds, R_Crds, T_Crds]
  #     x_Coords = [resname, resid, bb_crds, ca_crd, cg_crd, avg_crd]  

    # Check for missing residues
    for idx, R_Seq in enumerate(Input[5]):
      if R_Seq is None:
        print('  # R-spine Warning: Missing resid: '+Input[0]+' '+str(idx+1))
        return None
    for idx, C_Seq in enumerate(Input[6]):
      if C_Seq is None:
        print('  # C-spine Warning: Missing resid: '+Input[0]+' '+str(idx+1))
        return

  # R-spine calculation: use 2nd derivative to measure curvature of plot
  # each residue has 2 points, average of sidechain AND of Calpha. 
  # order of residues are given in the sequence input file
#    Refa = [np.asarray(R_Seq[ArrayCent(len(R_Seq))][3]) for R_Seq in self.ref_r]
#    Refc = [np.asarray(R_Seq[ArrayCent(len(R_Seq))][5]) for R_Seq in self.ref_r]

    Rsa  = [np.asarray(R_Seq[ArrayCent(len(R_Seq))][3]) for R_Seq in Input[5]]
    Rsc  = [np.asarray(R_Seq[ArrayCent(len(R_Seq))][5]) for R_Seq in Input[5]]
    Rsac = Rsa+Rsc
    Crd_r = np.asarray(zip(*Rsac))

    Rs_curv = []
    for x in range(3):
      m2, n1, c = np.polyfit(range(len(Crd_r[x])), Crd_r[x][:], 2, full=False)
      Rs_curv.append(m2*2)
    r_curv_m2 = VecMag(Rs_curv)

    r_fn2 = [LsqFit(range(len(Crd_r[x])), Crd_r[x], 2) for x in range(3)]
    R_Pts = [np.asarray([f(x) for f in r_fn2]) for x in range(len(Crd_r[x]))]

    r_curv = CalcCurvature(R_Pts)

  # C-spine calculation: use 2nd derivatives to measure curvature of plot
  # for calculation use 5 of the 8 C-spine residues defined in Roskoski 2016,
  # since the 3 of the residues are really not aligned linearly but more on
  # the side flanking the 5 residues stacking with ATP in 1ATP. The order of 
  # the residues are given in the sequence input file, although not important.
#    Ref = [np.asarray(T_Seq[ArrayCent(len(T_Seq))][5]) for T_Seq in self.ref_c]
    Css = [np.asarray(T_Seq[ArrayCent(len(T_Seq))][5]) for T_Seq in Input[6]]
    Csa = [np.asarray(T_Seq[ArrayCent(len(T_Seq))][3]) for T_Seq in Input[6]]
    Cs  = Csa+Css
    Crd_c = np.asarray(zip(*Cs))

    Cs_curv = []
    for x in range(3):
      m2,n1,c = np.polyfit(range(len(Crd_c[x])), Crd_c[x][:], 2, full=False)
      Cs_curv.append(m2*2)
    c_curv_m2 = VecMag(Cs_curv_8)

    c_fn2 = [LsqFit(range(len(Crd_c[x])), Crd_c[x], 2) for x in range(3)]
    C_Pts = [np.asarray([f(x) for f in c_fn2]) for x in range(len(Crd_c[x]))]

    c_curv = CalcCurvature(C_Pts)


    return [Input[0], r_curv, r_curv_m2, c_curv, c_curv_m2]
#    return [Input[0], r_curv, c_curv]


##########################################################################
# 3D-line curvature calculation. Input as array of 3D points
# only output the mid-point curvature (scalar) of the curvature vectors to
# represent the curvature
# result from the mid-point Acceleration (scalar) of the Acceleratio vectors
# is equivalent to the derivative of 2nd-order polyfit 1st slope parameter 
#   d(mx^2) = 2m (scalar)
# https://en.wikipedia.org/wiki/Curvature
# https://stackoverflow.com/questions/28269379/curve-curvature-in-numpy
def CalcCurvature( Array_3D ):

  a = Array_3D
  x, y, z = zip(*a)
  # calculate the derivatives of each variable and put back as Velocity
  dx_dt = np.gradient(x)
  dy_dt = np.gradient(y)
  dz_dt = np.gradient(z)
  velocity = np.array( [ [ dx_dt[i], dy_dt[i], dz_dt[i] ] 
                           for i in range(dx_dt.size) ] )

  # calculate Speed from Velocity
  ds_dt = np.sqrt( dx_dt**2 + dy_dt**2 + dz_dt**2 )

  # calculate tangent from Speed
  tangent = np.array( [1/ds_dt]*3 ).transpose() * velocity

  # calculate derivative of tangents
  tangent_x = tangent[:,0]
  tangent_y = tangent[:,1]
  tangent_z = tangent[:,2]
  d_tangt_x = np.gradient(tangent_x)
  d_tangt_y = np.gradient(tangent_y)
  d_tangt_z = np.gradient(tangent_z)
  
  dT_dt = np.array([[d_tangt_x[i],d_tangt_y[i],d_tangt_z[i]]
                     for i in range(d_tangt_x.size)])
  length_dT_dt = np.sqrt(d_tangt_x**2 + d_tangt_y**2 + d_tangt_z**2)
  normal = np.array([1/length_dT_dt]*3).transpose() * dT_dt

  # calculate 2nd derivatives for each component
  d2s_dt2 = np.gradient(ds_dt)
  d2x_dt2 = np.gradient(dx_dt)
  d2y_dt2 = np.gradient(dy_dt)
  d2z_dt2 = np.gradient(dz_dt)

  # calculate curvature
  a = (d2z_dt2*dy_dt - d2y_dt2*dz_dt)**2
  b = (d2x_dt2*dz_dt - d2z_dt2*dx_dt)**2
  c = (d2y_dt2*dx_dt - d2x_dt2*dy_dt)**2
  d = (dx_dt**2 + dy_dt**2 + dz_dt**2)**3
  curvature = np.sqrt( (a+b+c)/d )

  t_component = np.array([d2s_dt2]*3).transpose()
  n_component = np.array([curvature * ds_dt**2]*3).transpose()

  acceleration = t_component*tangent + n_component*normal

  # Take the curvature vector of middle of curve as representative
  if len(curvature) % 2 == 0:
    mid  = len(curvature)/2
    take = (curvature[mid-1] + curvature[mid])/2
  else:
    mid = len(curvature)-1/2
    take = curvature[mid]

  curv_mag = VecMag(take)
  return curv_mag


##########################################################################
#
def RCSpinesMeasure( Ref_Coords, Tgt_Coords, Data, output ):

  # Input_Coords = [pdb_name, H_Crds, N_Crds, C_Crds, G_Crds, R_Crds, T_Crds]
  #     x_Coords = [resname, resid, bb_crds, ca_crd, cg_crd, avg_crd]  

  # Get the R/C-spine Center
  pRC = RCSpines(ref_r=Ref_Coords[5], ref_c=Ref_Coords[6])
  Ref = pRC(Ref_Coords)

  # Create R/C-spine object for MPI
  mpi = multiprocessing.Pool(processes = multiprocessing.cpu_count())

#  Tmp = [pRC(Tgt) for Tgt in Tgt_Coords]
  Tmp = mpi.map(pRC, Tgt_Coords)
  mpi.close()
  mpi.join()

  # Tgt = [pdb_name, r_curv, r_curv_m2, c_curv, c_curv_m2]
  Tgt_List = [x for x in Tmp if x is not None]
  print('\n ## R-C Spine return: {0}\n'.format(len(Tgt_List)))
  CollectSpines(Ref, Tgt_List, Data)


##########################################################################

def CollectSpines(Ref, Tgt_List, Data):
  Tgt_List.append(Ref)
  for Tgt in Tgt_List:
    Data[Tgt[0]]['r_curv']    = Tgt[1]
    Data[Tgt[0]]['r_curv_m2'] = Tgt[2]
    Data[Tgt[0]]['c_curv']    = Tgt[3]
    Data[Tgt[0]]['c_curv_m2'] = Tgt[4]

#    Data[Tgt[0]]['c_sp_6']  = Tgt[5]
#    Data[Tgt[0]]['c_sp_8']  = Tgt[6]


##########################################################################
#
#   Peter Man-Un Ung @ MSSM
#
#   v1.0    17.03.10
#   v2.0    17.04.12    bugfix; deal with NoneType parameters
#   v3.0    17.08.19    change to curvature calculation for spines. Use COM
#                       of residue sidechain and backbone for R-/C-spines
#
