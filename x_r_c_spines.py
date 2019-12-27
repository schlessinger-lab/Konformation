#!/usr/bin/python

import re,os,glob,sys
import numpy as np
from pathos import multiprocessing
from tqdm import tqdm

from x_helix_axis  import *
from numpy.linalg import norm

##########################################################################

  ### For R-/C-spine curvature calculation, we assume the stacking of residues
  ### is a simple 2nd-order curve and not anything higher order (3rd, 4th, etc)
  ### and we measure how linear the stacking of the residues by calculating 
  ### the curvature of a simple 2nd-order curve in a 3-dimensional (Cartesian) 
  ### space as a scalar number for comparison
  ### In the spine curvature calculation, use the side chain and backbone,
  ### generate the 2nd-order curve function and populate the points by
  ### using more number

##########################################################################
class RCSpines( object ):
  
  ## use of reference R-/C-spine data is depricated
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
        print('\n  #2# R-spine Warning: Missing resid: '+Input[0]+' '+str(idx+1))
        return None
    for idx, C_Seq in enumerate(Input[6]):
      if C_Seq is None:
        print('\n  #2# C-spine Warning: Missing resid: '+Input[0]+' '+str(idx+1))
        return None

  # R-spine calculation: 
  # order of residues are given in the sequence input file
    # get the (x,y,z) polynominal functions and vectors (points) of a curve
    # hard-code to get 21 points from formula
    # same result with 21 points vs. np.arange(0., len(Crd_r[0]), 0.25)
    Rbb  = [np.asarray(R_Seq[ArrayCent(len(R_Seq))][3]) for R_Seq in Input[5]]
    Rsc  = [np.asarray(R_Seq[ArrayCent(len(R_Seq))][5]) for R_Seq in Input[5]]
    Rsac = Rsc+Rbb

    Crd_r = np.asarray(zip(*Rsac))

    r_fn2 = [LsqFit(range(len(Crd_r[x])), Crd_r[x], 2) for x in range(3)]
    R_Pts = [np.asarray([f(x) for f in r_fn2]) for x in range(21)]

    # CalCurvature 1, 2, 3 have different codes but same equation, same result
    r_curv1 = CalcCurvature1(R_Pts)	# Generalized form
#    r_curv2 = CalcCurvature2(R_Pts)	# Parametric Cartesian form
#    r_curv3 = CalcCurvature3(R_Pts)	# Long Cartesian form
#    print('{0:8.5f} {1:8.5f} {2:8.5f}'.format(r_curv1, r_curv2, r_curv3))

  # C-spine calculation: 
  # for calculation use 5 of the 8 C-spine residues defined in Roskoski 2016,
  # since the 3 of the residues are really not aligned linearly but more on
  # the side flanking the 5 residues stacking with ATP in 1ATP. The order of 
  # the residues are given in the sequence input file, although not important.
    Cbb = [np.asarray(T_Seq[ArrayCent(len(T_Seq))][3]) for T_Seq in Input[6]]
    Csc = [np.asarray(T_Seq[ArrayCent(len(T_Seq))][5]) for T_Seq in Input[6]]
    Cs  = Csc+Cbb

    Crd_c = np.asarray(zip(*Cs))

    c_fn2 = [LsqFit(range(len(Crd_c[x])), Crd_c[x], 2) for x in range(3)]
    C_Pts = [np.asarray([f(x) for f in c_fn2]) for x in range(21)]

    c_curv1 = CalcCurvature1(C_Pts)
#    c_curv2 = CalcCurvature2(C_Pts)

    return [Input[0], r_curv1, c_curv1]


##########################################################################
## Curvature of a 3D curve based on curvature formula in Wikipedia
## https://en.wikipedia.org/wiki/Curvature
## no need for the complex tangent/acceleration calculation in
## https://stackoverflow.com/questions/28269379/curve-curvature-in-numpy
#def CalcCurvature2( Curve ):
  
  ## Generalized curvature expression, independent of the nth-dimension
  ## first calculate 1st and 2nd derivatives of list of curve points
  ## double vertical bar = norm (length) of a vector (not unit vector)
  ## curvature k = || r' x r'' || / || r' ||**3
#  dr_dt     = np.gradient( Curve, axis=0 )
#  d2r_dt2   = np.gradient( dr_dt, axis=0 )
#  curvature = [ norm(np.cross(dr, d2r)) / (norm(dr)**3) 
#                for dr, d2r in zip(dr_dt, d2r_dt2) ]

  # Take the curvature vector of middle of curve as representative
#  if len(curvature) % 2 == 0:
#    mid  = len(curvature)/2
#    take = (curvature[mid-1] + curvature[mid])/2
#  else:
#    mid = (len(curvature)-1)/2
#    take = curvature[mid]

#  return curvature[mid]


#####################
## Curvature calculation through Parametric 3-D (only) Cartesian coodinates
def CalcCurvature3( Curve ):
  # older way to calculate 1st derivatives of list of curve pts in 3-dimension
  x, y, z = zip(*Curve)    # dx_dt = np.gradient(x)
  dr_dt = np.asarray( zip(np.gradient(x), np.gradient(y), np.gradient(z)) )
  d2r_dt2   = np.gradient( dr_dt, axis=0 )

  ## Curvature of parametric Cartesian coordinates by r(i) = (x(i),y(i),z(i))
  ## curvature k = sqrt( (z"y' - y"z')^2 + (x"z' - z"x')^2 + (y"x' - x"y')^2 )/
  ##               sqrt( (x'^2 + y'^ + z'^2)**3 )
  dr, d2r = dr_dt, d2r_dt2
  curvature = []
  for i in range(np.size(dr_dt, axis=0)):
    x = (( d2r[i][2]*dr[i][1] - d2r[i][1]*dr[i][2] )**2 +
         ( d2r[i][0]*dr[i][2] - d2r[i][2]*dr[i][0] )**2 +
         ( d2r[i][1]*dr[i][0] - d2r[i][0]*dr[i][1] )**2 )
    y = ( dr[i][0]**2 + dr[i][1]**2 + dr[i][2]**2 )**3
    curvature.append( np.sqrt(x/y) )

  # Take the curvature vector of middle of curve as representative
  if len(curvature) % 2 == 0:
    mid  = len(curvature)/2
    take = (curvature[mid-1] + curvature[mid])/2
  else:
    mid = (len(curvature)-1)/2
    take = curvature[mid]

  return curvature[mid]


##########################################################################
# 3D-line curvature calculation. Input as array of 3D points
# only output the mid-point curvature (scalar) of the curvature vectors to
# represent the curvature
# result from the mid-point Acceleration (scalar) of the Acceleration vectors
# is equivalent to the derivative of 2nd-order polyfit 1st slope parameter 
#   d(mx^2) = 2m (scalar)	which is also the tangent of the curve
# https://en.wikipedia.org/wiki/Curvature
# https://stackoverflow.com/questions/28269379/curve-curvature-in-numpy
def CalcCurvature1( Array_3D ):

  a = Array_3D
  x, y, z = zip(*a)
  # calculate the derivatives of each variable and put back as Velocity
  dx_dt = np.gradient(x)
  dy_dt = np.gradient(y)
  dz_dt = np.gradient(z)
  velocity = np.array( [ [ dx_dt[i], dy_dt[i], dz_dt[i] ] 
                           for i in range(dx_dt.size) ] )

  # calculate 2nd derivatives for each component
  d2x_dt2 = np.gradient(dx_dt)
  d2y_dt2 = np.gradient(dy_dt)
  d2z_dt2 = np.gradient(dz_dt)
  d2r_dt2 = np.array( [ [ d2x_dt2[i], d2y_dt2[i], d2z_dt2[i] ]
                          for i in range(d2x_dt2.size) ] )

  # calculate curvature
  a = (d2z_dt2*dy_dt - d2y_dt2*dz_dt)**2
  b = (d2x_dt2*dz_dt - d2z_dt2*dx_dt)**2
  c = (d2y_dt2*dx_dt - d2x_dt2*dy_dt)**2
  d = (dx_dt**2 + dy_dt**2 + dz_dt**2)**3
  curvature = np.sqrt( (a+b+c)/d )


  # calculate Speed from Velocity
#  ds_dt   = np.sqrt( dx_dt**2 + dy_dt**2 + dz_dt**2 )
#  d2s_dt2 = np.gradient(ds_dt)

  # calculate tangent from Speed
#  tangent = np.array( [1/ds_dt]*3 ).transpose() * velocity

  # calculate derivative of tangents
#  tangent_x = tangent[:,0]
#  tangent_y = tangent[:,1]
#  tangent_z = tangent[:,2]
#  d_tangt_x = np.gradient(tangent_x)
#  d_tangt_y = np.gradient(tangent_y)
#  d_tangt_z = np.gradient(tangent_z)
  
#  dT_dt = np.array([[d_tangt_x[i],d_tangt_y[i],d_tangt_z[i]]
#                     for i in range(d_tangt_x.size)])
#  length_dT_dt = np.sqrt(d_tangt_x**2 + d_tangt_y**2 + d_tangt_z**2)
#  normal = np.array([1/length_dT_dt]*3).transpose() * dT_dt

#  t_component = np.array([d2s_dt2]*3).transpose()
#  n_component = np.array([curvature * ds_dt**2]*3).transpose()
#  acceleration = t_component*tangent + n_component*normal

  # Take the curvature vector of middle of curve as representative
  if len(curvature) % 2 == 0:
    mid  = len(curvature)/2
    take = (curvature[mid-1] + curvature[mid])/2
  else:
    mid = (len(curvature)-1)/2
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
  Tmp = [x for x in tqdm(mpi.imap_unordered(pRC, Tgt_Coords),total=len(Tgt_Coords))]
  mpi.close()
  mpi.join()

  # Tgt = [ pdb_name, r_curv, c_curv ]
  Tgt_List = [x for x in Tmp if x is not None]
  print('\n ## R-C Spine return: {0}\n'.format(len(Tgt_List)))
  CollectSpines(Ref, Tgt_List, Data)


##########################################################################

def CollectSpines(Ref, Tgt_List, Data):
  Tgt_List.append(Ref)
  for Tgt in Tgt_List:
    Data[Tgt[0]]['r_curv'] = Tgt[1]
    Data[Tgt[0]]['c_curv'] = Tgt[2]

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
#   v4.0    17.08.29    change curvature calculation to simpler format
#
