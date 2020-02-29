#!/usr/bin/env python3

import re,sys
import numpy as np

from tqdm import tqdm
from numpy.linalg import norm
from pathos import multiprocessing
np.seterr(invalid='ignore')

from aa_residue import AA

##########################################################################
## Compare the Helix axis and distance between the Reference and Target axes
## Run using MPI then print out the results into an output file
def HelixMeasurements( Ref_Coords, Tgt_Coords, Data, parm, output ):
  # Input_Coords = [pdb_name, H_Crds, N_Crds, C_Crds, G_Crds, R_Crds, T_Crds]
  #     x_Coords = [resname, resid, bb_crds, ca_crd, cg_crd, avg_crd, cb_crd] 

  print('##################################################################\n')

  # Create helix object for MPI
  Ref = HelixAxis(Ref_Coords)

  if parm['MPICPU'][0] == 1:
    Tmp = [ HelixAxis(Tgt) for Tgt in Tgt_Coords ]
  else:
    if parm['MPICPU'][0] == 0:
      mpi_cpu = multiprocessing.cpu_count()
    else:
      mpi_cpu = parm['MPICPU'][0]
    mpi   = multiprocessing.Pool(mpi_cpu)
    Tmp = [x for x in tqdm(mpi.imap_unordered(HelixAxis, Tgt_Coords), total=len(Tgt_Coords))]
    mpi.close()
    mpi.join()

  # Tgt = [ pdb_id, res_id, axis, cg_nom, cg_vec, sc_vec,
  #         sc_pres, cg_pres, curve, phi, psi, r_median, r_std,
  #         Reg2 ]

  ## Have to leave None in Reg2 to keep Reg2 to have same number as PDB_Coords
  ## for DomainDistances(), which needs Reg2
#  Tgt_List = [x for x in Tmp if x is not None]
  Tgt_List = Tmp
  print('\n ## Helix Axis return: {0}\n'.format(len(Tgt_List)))
  CollectHelix(Ref, Tgt_List, Data)

  # Extract the 2nd-order regression data for other use
  # last item in Tgt_List is the Reference item.. pop
  RefReg2 = Ref[-1]
  Reg2    = [ Tgt[-1] for Tgt in Tgt_List[:-1] ]
  return RefReg2, Reg2


#########################################################################
def CollectHelix( Ref, Tgt_List, Data ):
  # { [name, resid, vec, cg_nom, cg_vec, sc_vec, sc_pres, cg_pres, 
  #                 curve, phi, psi, r_median, r_std, Reg2] }

  Tgt_List.append(Ref)
  for Tgt in Tgt_List:
    if Tgt[0] is not None:
      Data[Tgt[0]]['h_sc_x'] = Tgt[6]
      Data[Tgt[0]]['h_cg_x'] = Tgt[7]
      Data[Tgt[0]]['h_curv'] = Tgt[8]
      Data[Tgt[0]]['h_phi']  = Tgt[9]
      Data[Tgt[0]]['h_psi']  = Tgt[10]
      Data[Tgt[0]]['r_medi'] = Tgt[11]
      Data[Tgt[0]]['r_std']  = Tgt[12]

      if Tgt[2] is not None:
        Data[Tgt[0]]['h_axvc'] = np.dot(Tgt[2]/VecMag(Tgt[2]),
                                        Ref[2]/VecMag(Ref[2]))
      if Tgt[3] is not None:
        Data[Tgt[0]]['h_norm'] = np.dot(Tgt[3]/VecMag(Tgt[3]),
                                        Ref[3]/VecMag(Ref[3]))
      if Tgt[4] is not None:
        Data[Tgt[0]]['h_cgvc'] = np.dot(Tgt[4]/VecMag(Tgt[4]),
                                        Ref[4]/VecMag(Ref[4]))
      if Tgt[5] is not None:
        Data[Tgt[0]]['h_scvc'] = np.dot(Tgt[5]/VecMag(Tgt[5]),
                                        Ref[5]/VecMag(Ref[5]))


##########################################################################
## Calculate all helix properties, target residue at center of sequence
def HelixAxis( Input ):
  # Input_Coords = [pdb_name, H_Crds, N_Crds, C_Crds, G_Crds, R_Crds, T_Crds]
  #       x_Crds = [resname, resid, bb_crds, ca_crd, cg_crd, avg_crd, cb_crd]
  None_Coord = [None, None, None, None, None, None, None, None, None,
                None, None, None, None ]
  
  pdb_id, Pre_Coords = Input[0:2]

  # Check for missing residues
  if Pre_Coords is None:
    print('  # Helix Warning: Too short to calculate: \033[31m{0}]033[0m'.format(pdb_id))
    return None_Coord
  for idx, Seq in enumerate(Pre_Coords):
    if Seq is None:
      print('\n  #2# Helix Warning: Missing residue: {0}\t{1}'.format(pdb_id,idx+1))
      return None_Coord


  Center = Pre_Coords[ArrayCent(len(Pre_Coords))]
  res_id = AA(Center[0])+str(Center[1])
  
  # reformat the backbone atoms N,CA,C array (bb_coords)
  try:
    Coords = np.asarray(sum(list(zip(*Pre_Coords))[2],[]))
  except TypeError:
    print('\n  #2# Helix Warning: Cannot read data - TypeError: '+pdb_id)
    return None_Coord

  if len(Coords) < 7:
    print('\n  #2# Helix Warning: Helix too short to calculate: \033[31m{0}\033[0m'.format(pdb_id))
    return None_Coord

  # Calculate helix normal vector and helix center line (1st/2nd-regression),
  # helix radii over the length using 2nd-order axis, and spherical angles 
  # of 1st-order helix axis
  axis, curve, Reg1, Reg2 = CalculateHelixAxis(Coords)
  phi, psi = SphericalAngles(axis)
  r_median, r_std = HelixRadius(Coords, Reg2, pdb_id)

  # Calculate CB, CG, SC vectors; C-axis/Glu-cg normal vector (cross product)
  # SC vector will be available as long as CB is present
  # sc_pres = h_sc_x = cg_vec.sc_vec	cg_pres = h_cg_x = cg_vec.cb_vec
  cg_vec, cg_nom, sc_vec, sc_pres, cg_pres = None, None, None, None, None
  if Center[4] is not None:
    cg_vec  = np.array( Center[4] - Reg2[int(ArrayCent(len(Coords)))] )
    sc_vec  = np.array( Center[5] - Reg2[int(ArrayCent(len(Coords)))] )
    cb_vec  = np.array( Center[6] - Reg2[int(ArrayCent(len(Coords)))] )

    cg_nom  = np.cross(axis/VecMag(axis), cg_vec/VecMag(cg_vec))
    sc_pres = np.dot(cg_vec/VecMag(cg_vec), sc_vec/VecMag(sc_vec))
    cg_pres = np.dot(cg_vec/VecMag(cg_vec), cb_vec/VecMag(cb_vec))
  else:
    # If CB is available, substitute CG with CB for CG-vector
    if Center[6] is not None:
      print('\n  #1# Helix Warning: No CG, use "CB" for Ax-Cg vector: '+pdb_id)
      cb_vec = np.array( Center[6] - Reg2[int(ArrayCent(len(Coords)))] )
      sc_vec = np.array( Center[5] - Reg2[int(ArrayCent(len(Coords)))] )
      cg_vec = cb_vec

      cg_nom  = np.cross(axis/VecMag(axis), cg_vec/VecMag(cg_vec))
      if Center[5] is not None:
        sc_pres = np.dot(cg_vec/VecMag(cg_vec), sc_vec/VecMag(sc_vec))
    else:
      print('\n  #1# Helix Warning: No CG for Ax-Cg vector: '+pdb_id)


  if curve is None:
    print('\n  #1# Helix Warning: Helix residues < 5, skip "curv": '+pdb_id)

  return  [ pdb_id, res_id, axis, cg_nom, cg_vec, sc_vec, 
            sc_pres, cg_pres, curve, phi, psi, r_median, r_std, 
            Reg2 ]


#########################################################################
## Calculate the helix axis using coordinates supplied, calculate 1st- and 2nd-
## order regression curves to represent helix axis. Calculate the helix 
## curvature centering at conserved Glu
def CalculateHelixAxis( Input ):

  Coords = np.asarray(Input)
  count  = len(Coords)

  # Linear regression on Cartesian Coordinates to estimate helix axis coords.
  # Use moving sets of points on both end to average out regression error
  # iterate range(3) to calculate x,y,z coordinates separately.
  # Use minimium of 7 residues (> 21 points),
##  # If available residue is < 5 (< 15 points), helix is incomplete and calc
##  # will be inaccurate, skip h_vec and h_curve with None
  if count >= 15:
    posit  = 6                # 3 atoms (N,C,CA) = 1 residue
  else:
    posit  = 1
  xcount = count - posit    # reduced number of points to do LSQ

  ## maximum of 10 tries to calculate Reg2 until it yield result
  ## this is to counter a bug with numpy ployfit where it fails randomly
  ## when running a lot of things but doing fine when only running a few things
  for idx in range(10):
    Fn1Pts, Fn2Pts = [], []
    for m in range(0,posit):
      Fn1 = [LsqFit(list(range(xcount)), Coords[m:m-posit, x],1) for x in range(3)]
      Fn1Pts.append( [np.asarray([f(x) for f in Fn1]) for x in range(count) ])
      Fn2 = [LsqFit(list(range(xcount)), Coords[m:m-posit, x],2) for x in range(3)]
      Fn2Pts.append( [np.asarray([f(x) for f in Fn2]) for x in range(count) ])

    H_Curv = [ CalcCurvature2(Fn2) for Fn2 in Fn2Pts ]
#    h_curv = ( H_Curv[2]+H_Curv[3] )/2
#    hc_std = np.std(H_Curv[2:4])

    Reg1 = np.mean(Fn1Pts, axis=0)
    Reg2 = np.mean(Fn2Pts, axis=0)
  Start, End, Center = Reg1[0], Reg1[-1], Reg2[ ArrayCent(count) ]

  if posit > 1:   
    vec   = (End-Start)/VecMag(End-Start)
#    curve = VectorAngle( (Center-Start), (End-Center) )
    curve = ( H_Curv[2]+H_Curv[3] )/2
  else:
    vec   = Center/VecMag(Center)
    curve = None

  return [ vec, curve, Reg1, Reg2 ]


##########################################################################
## Calculate the spherical angles of a vector relative to z-axis
## Phi as angle (y-axis); Psi as dihedral angle (z-axis) between 2 vectors
def SphericalAngles(vec):

  # phi for vector to xz-plane/x_axis angle
  x_axis = [1,0,0] 
  norm   = VecMag([vec[0],vec[1],0])
  phi    = np.arccos(np.dot([vec[0],vec[1],0], x_axis)/norm)*180/np.pi

  # psi for vector to yz-plane/z-axis angle
  z_axis = [0,0,1]
  h_norm = VecMag(vec)
  psi    = np.arccos(np.dot(vec, z_axis)/h_norm)*180/np.pi

  return phi, psi


##########################################################################
## Calculate average radius of helix using 2nd-order regression as helix center
def HelixRadius(Coords, Reg2Pts, infile):

  count = len(Coords)
  if count != len(Reg2Pts):
    sys.exit('\n  \033[31m#2# Helix FATAL:\033[0m No. of regression points does not match number of coord points: \033[35m{0}\033[0m'.format(Coords[0]))

  # Cylindrical coodinates of helix
  Dist  = [VecMag(Coords[i]-Reg2Pts[i]) for i in range(count)]
  d     = np.poly1d(np.polyfit(list(range(count)), Dist, 1, full=False))
  median, stdev = np.median(Dist), np.std(Dist)

  return median, stdev


##########################################################################
## Curvature of a 3D curve based on curvature formula in Wikipedia
## https://en.wikipedia.org/wiki/Curvature
## no need for the complex tangent/acceleration calculation in
## https://stackoverflow.com/questions/28269379/curve-curvature-in-numpy
def CalcCurvature2( Curve ):

  ## Generalized curvature expression, independent of the nth-dimension
  ## first calculate 1st and 2nd derivatives of list of curve points
  ## double vertical bar = norm (length) of a vector (not unit vector)
  ## curvature k = || r' x r'' || / || r' ||**3
  dr_dt     = np.gradient( Curve, axis=0 )
  d2r_dt2   = np.gradient( dr_dt, axis=0 )
  curvature = [ norm(np.cross(dr, d2r)) / (norm(dr)**3) 
                for dr, d2r in list(zip(dr_dt, d2r_dt2)) ]

  # Take the curvature vector of middle of curve as representative
  if len(curvature) % 2 == 0:
    mid  = int( len(curvature)/2 )
    take = (curvature[mid-1] + curvature[mid])/2
  else:
    mid = int( (len(curvature)-1)/2 )
    take = curvature[mid]

  return take


##########################################################################
## n-order Least-square fit of 2 arrays, x and y, return an object/def as
## a nst-order polynomial function f(x) = mx+c, or f(x) = mx^2+nx+c, etc
## i.e.                            f(1) = 15       f(5) = 25
def LsqFit( X, Y, order=1 ):
  return np.poly1d( np.polyfit( X, Y, order, full=False))


##########################################################################
# Find the center element in an array. The result is 1 less than correct
# answer -- Python numbering starts from 0. i.e. 9-element array is [0,...,8]
# center is 4, not 5
def ArrayCent( count ):
  if count % 2 == 0:
    center = count/2-1
  else:              
    center = ((count-1)/2)
  return int(center)


##########################################################################
## Magnitude of a vector, equivalent to np.linalg.norm( v )
## this algorithm is much faster than np.linalg.norm function
def VecMag( v ):
  mag = np.sqrt((np.array(v)**2).sum(-1))
  return mag

## Angle between vectors
def VectorAngle( v1, v2 ):
  ang = np.arccos( np.dot(v1,v2)/(VecMag(v1)*VecMag(v2)) ) * 180/np.pi
  return ang

## Distance between 2 points
def Distance( v1, v2 ):
  return VecMag( v2 - v1 )

## Distance between vectors
def VectorDistance(start1, v1, start2, v2):
  cross_prod = np.cross(v1,v2)
  mx         = VecMag(cross_prod)
  norm       = cross_prod/mx
  diff       = start1 - start2
  dist       = np.fabs( np.dot(norm,diff) )
  return dist


##########################################################################
##
##  Peter M.U. Ung  @   MSSM
##
##  v0.1    17.01.08
##  v0.2    17.01.29
##  v0.3    17.01.30
##
##  v1.0    17.02.02
##  v1.1    17.02.05
##  v2.0    17.02.21	- add data collection object
##  v2.1    17.04.12    bug fix; check for missing residues
##  v3.0    17.04.28    add H-axis/Glu-cg normal vector comparison
##  v4.0    17.04.28    fixed a big error in Psi spherical angle calculation
##  v5.0    17.05.08    add CB and SC vectors; add presence helix SC
##  v6.0    17.08.29    update helix curvature calculation method
##  v7.0    17.10.09    cleanup and update list to keep same list_length
##  v8.0    18.03.12    bug fix for helix algorithm
##
##########################################################################
