#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 18:52:26 2019

@author: pmung
"""

import sys,os
import re
#import rpy2
import time
import pickle
import numpy as np
import pandas as pd
import mdtraj as md
from tqdm import tqdm
from pathos import multiprocessing
from argparse import ArgumentParser
np.seterr(invalid='ignore')

lib_dir = '~/Dropbox (Schlessinger lab)/9_scripts/3_program/structures/4_Konformation/z_database/'
wrk_dir = '~/Dropbox (Schlessinger lab)/z_others/8_strada/'

lib_dir = '~/Dropbox/9_scripts/3_program/structures/4_Konformation/z_database/'
wrk_dir = '~/Dropbox/z_others/8_strada/'
##########################################################################
def main():
  args = UserInput()

  args.tmpl_file = wrk_dir+'2_md/strada_cido.prot.1atp.pdb'
  args.traj_file = wrk_dir+'2_md/strada_cido.1.200ps.dcd'
  args.outpref   = 'test'
  args.b3k = 39
  args.dfg = 152
  args.c_glu = 57

  ## reference structure must start at resid 1. Modified ref is hardcoded here
  if not os.path.isfile(lib_dir+'1ATP.mdtraj.pdb'):
    sys.exit('\n    ERROR: Reference structure "1ATP.mdtraj.pdb" not found\n')
  else:
    ref_file = lib_dir+'1ATP.mdtraj.pdb'
    ref_pkl  = lib_dir+'1ATP.mdtraj.pkl'
    ref_dfg  = 171
    ref_b3k  = 59
    ref_c_glu= 78

######################

  ## get reference PDB structure 1ATP.pdb coordinates dataframe
  print('# Reading in reference file: '+ref_file)
  if not ref_pkl or not os.path.isfile(ref_pkl):
    ref    = md.load_pdb(ref_file)
    ref_cd = ExtractCoords(dfg=ref_dfg, b3k=ref_b3k, c_glu=ref_c_glu, pkl=ref_pkl)
    ref_df = CalculateMetrics( ref_cd(ref) )
  ## skip calculation if data is already stored in pickle
  else:
    print('  ## INFO: Read structural residue coords from: {0}\n'.format(ref_pkl))
    ref_df = CalculateMetrics( pickle.load( open(ref_pkl, 'rb') ) )

######################
  ## load trajectory file(s) with MDtraj, can be multiple traj files at once
  traj = []
  print('# Reading in trajectory file(s)...')
  start = time.perf_counter()
  if not args.pkl or not os.path.isfile(args.pkl):
    start2= time.perf_counter()
    TrjIn = ReadTraj(top=args.tmpl_file) 
    if re.search(r'dcd$|nc$|crd$|xtc$', args.traj_file):
      traj = TrjIn(args.traj_file)
    else:
      traj_list = filter(None, (l.rstrip() for l in open(args.traj_file, 'r')
                                if l is not re.search(r'^#', l)))
      mpi  = multiprocessing.Pool(processes = multiprocessing.cpu_count())
      traj = md.join( mpi.imap(TrjIn, traj_list) )
      mpi.close()
      mpi.join()
    end2 = time.perf_counter()
    print('  ## Time to load trajectory: {0:.1f} ms for {1} frames\n'.format(
             (end2-start2)*1000, len(traj)) )

  ## superpose all frames to template structure pre-superposed to ref 1ATP.pdb
    if args.superp:
      print('# Applying superposition to trajectory with: '+args.superp)
      tmpl = md.load_pdb(args.tmpl_file)
      traj = traj.superpose(tmpl, atom_indices=args.superp, parallel=True)

    ## get trajectory coordinates dataframe
    print('# Extracting structural matrics from trajectory...')
    start  = time.perf_counter()
    trj_cd = ExtractCoords(dfg=args.dfg, b3k=args.b3k, c_glu=args.c_glu, pkl=args.pkl)
    trj_df = CalculateMetrics( trj_cd(traj) )
             
  ## skip calculation if data is already stored in pickle
  else:
    print('  ## INFO: Read structural residue coords from: {0}\n'.format(args.pkl))
    trj_df = CalculateMetrics( pickle.load( open(args.pkl, 'rb') ) )

  end    = time.perf_counter()
  print('## Total time to get traj descriptors: {0:.1f} ms for {1} frames'.format(
        (end-start)*1000, len(trj_df)))
  del traj    # save memory
  print('\n#########################################\n')

######################
######################
  ## calculate structural metrics from coordinates, then print out raw output
  print('# Calculating structural matrics from coordinates...')
  start  = time.perf_counter()
  mat_df = CompareMetrics(trj_df, ref_df)
  end    = time.perf_counter()
  print('## Total time to compare descriptors: {0:.1f} s for {1} frames'.format(
        (end-start)*1000000, len(mat_df)))
  print('\n#########################################\n')

  mat_df.to_csv(args.outpref+'.csv', sep=',')
  mat_df.to_hdf(args.outpref+'.h5', key='df', mode='w')

#####################
  ## use Kinformation Random Forest Classifier to assign conformation/confidence
#  mat_df['conform'] = Kinformation(mat_df)
  mat_df.index.name = '#frame'

  print('\n#########################################\n')



##########################################################################
##########################################################################
## determine kinase conformation of array of frames, save "conf" and confidence
def Kinformation( mat_df ):
  Matr = ['conf', 'cidi',      'cido',      'codi',      'codo',      'wcd',
                  'cidi-conf', 'cido-conf', 'codi-conf', 'codo-conf', 'wcd-conf']
#  mat_df = ['p1p1x', 'p2p2x', 'v3v3x', 'dfg_st', 'h_cgvc',
#            'ang_NHs', 'ang_CHs', 'dist_NC', 'dist_NH', 'dist_CH']

  Cols = ['conf']
  f_df = pd.DataFrame(index=range(0, len(mat_df), columns=Matr))

  import rpy2
  library(randomForest)
  library(clusterSim)
  load('./chelixmod.rda')
  load('./dfgmod.rda')
  ## matrics dataframe
#pdb_id(2), Group(4), p1p1x(5), p2p2x(6), h_cgvc(7),
#ang_NHs(8), ang_CHs(9), dist_NC(10), dist_NH(11), dist_CH(12),
#2,4,5,6,7,8,10,11,13




  return f_df


##########################################################################
## calculate structural vectors, and recycle the angular/distance metrics
def CompareMetrics( trj_df, ref_df_orig ):
  Cols = ['p1p1x', 'p2p2x', 'v3v3x', 'dfg_st', 'h_cgvc',
          'ang_NHs', 'ang_CHs', 'dist_NC', 'dist_NH', 'dist_CH']
  c_df = pd.DataFrame(index=range(len(trj_df)), columns=Cols)

  # for vectorization, make ref_df and trj_df same dimension
  ref_df = pd.DataFrame(np.repeat(ref_df_orig.values, len(trj_df), axis=0))
  ref_df.columns = ref_df_orig.columns

  # transpose the input data for vectorization, but no need to transpose back
  # when exit - already in 1-D array format
  c_df.p1p1x   = VecDot([ np.array(list(trj_df['p1'])).T, 
                          np.array(list(ref_df['p1'])).T ])
  c_df.p2p2x   = VecDot([ np.array(list(trj_df['p2'])).T, 
                          np.array(list(ref_df['p2'])).T ])
  c_df.v3v3x   = VecDot([ np.array(list(trj_df['v3'])).T, 
                          np.array(list(ref_df['v3'])).T ])
  c_df.h_cgvc  = VecDot([ np.array(list(trj_df['cg_vec'])).T, 
                          np.array(list(ref_df['cg_vec'])).T ])
  c_df.dfg_st  = DFGState( c_df )

  c_df.ang_NHs = trj_df.ang_NHs
  c_df.ang_CHs = trj_df.ang_CHs
  c_df.dist_NC = trj_df.dist_NC
  c_df.dist_NH = trj_df.dist_NH
  c_df.dist_CH = trj_df.dist_CH

  return c_df

###############################################################################
def DFGState( df ):
  ## Model PDB has same DFG- config as template DFG-in:     'in'
  ## Model PDB has opposite DFG- config as template DFG-in: 'out'
  ## Model PDB has undefined DFG- config:                   'random'

  ## predefine a table of value for Dataframe
  ## numpy vectorized method, faster by 125-1000x
  def _conditions2( x, y ):
    dfg_in  = (x > 0.005) & (y > 0.050)
    dfg_out = (x <-0.125) & (y <-0.125)

    dfg = ['weird']*len(x)
    for i in range(len(x)):
      if dfg_in[i]:
        dfg[i] = 'in'
      if dfg_out[i]:
        dfg[i] = 'out'
    return dfg

  ## original conditions method
  def _conditions( df ):
    if df.p1p1x is None or df.p2p2x is None:
      return 'missing DFG'
    elif df.p1p1x > 0.005 and df.p2p2x > 0.05:
      return 'in'
    elif df.p1p1x < -0.125 and df.p2p2x < -0.125:
      return 'out'
    else:
      return 'weird'

  df.dfg_st = _conditions2(df.p1p1x.to_numpy(), df.p2p2x.to_numpy())
#  df.dfg_st = df.apply(_conditions, axis=1)    # slower by 100x
  return df.dfg_st


##########################################################################
## Calculate the structural metrics using the list of coordinates
## Instead of using MPI looping over the coordinates, here use Dataframe + Numpy
## to 'Vectorize' the simple calculations (add, subtract, multiply, divide) to
## gain ~ 1,000-1,250x speed boost in comparsion to single processor operation.
## This vectorization can only be done to simple maths and cannot be done over
## complicated ones, e.g. linear regression, linear alegbra. For Numpy operations
## such as dot and cross products, have to break down the calculation back to they
## basic (+|-|x|/) operations.
def CalculateMetrics( Crd ):
  frames = len(list(Crd.coord_b3k.ca))
  ## predefine coord dataframe index/column (frame/atom) to maximize efficiency
  Cols = [ 'p1', 'p2', 'v3', 'cg_vec', 'ang_NHs', 'ang_CHs',
           'dist_NC', 'dist_NH', 'dist_CH', 'temp' ]
  m_df = pd.DataFrame(index=range(frames), columns=Cols)

  # to successfully put object of 'list of np.arrays' into rows of Dataframe cell,
  # convert the object into 'list' first, then into 'dictionary', then Dataframe
  d_df = pd.DataFrame({ 
    'dfg_d_ca': list(Crd.coord_dfg_d.ca), 'dfg_d_cb': list(Crd.coord_dfg_d.cb), 
    'dfg_d_cg': list(Crd.coord_dfg_d.cg), 'dfg_f_ca': list(Crd.coord_dfg_f.ca), 
    'dfg_f_cb': list(Crd.coord_dfg_f.cb), 'dfg_f_cg': list(Crd.coord_dfg_f.cg),
    'b3k_ca':   list(Crd.coord_b3k.ca),   'b3k_cb':   list(Crd.coord_b3k.cb), 
    'b3k_cg':   list(Crd.coord_b3k.cg),
    'c_glu_ca': list(Crd.coord_c_glu.ca), 'c_glu_cb': list(Crd.coord_c_glu.cb), 
    'c_glu_cg': list(Crd.coord_c_glu.cg), 'hlx_cent': list(Crd.coord_hlx.cn),
    }, index=range(frames) )
 
  print('#########################################\n')

####################
  ## Using Dataframe vectorization is ~ 250x faster than MPIx4, ~ 1000x looping
  print('# Calculate C-Glu vector (cg_vec)...')
  start = time.perf_counter()
  m_df.temp = VecGen( [d_df.c_glu_cg.to_numpy(), d_df.hlx_cent.to_numpy()] )
  m_df.cg_vec = m_df['temp'].div( VecMag( list(m_df['temp']) ) )
  end   = time.perf_counter()
  print(' {0:.1f} ms for {1} frames\n'.format((end-start)*1000, frames))
  del m_df['temp']

##################
  print('# Calculate DFG vectors (p1, p2, v3)...')
  start = time.perf_counter()
  Tmp = CalculateDFGVectors( 
                [d_df.dfg_d_cg.to_numpy(), d_df.dfg_d_ca.to_numpy(), 
                 d_df.dfg_f_ca.to_numpy(), d_df.dfg_f_cg.to_numpy()] )
  m_df.p1 = Tmp.p1
  m_df.p2 = Tmp.p2
  m_df.v3 = Tmp.v3
  end   = time.perf_counter()
  print(' {0:.1f} ms for {1} frames\n'.format((end-start)*1000, frames))

#####################
  print('# Calculate N-domain/C-helix angle (ang_NHs)...')
  start = time.perf_counter()
  m_df.ang_NHs = VectorAngle([ d_df.hlx_cent.to_numpy(), d_df.b3k_ca.to_numpy(),
                               d_df.hlx_cent.to_numpy(), d_df.c_glu_cg.to_numpy() ])
  end = time.perf_counter()
  print(' {0:.1f} ms for {1} frames\n'.format((end-start)*1000, frames))

#####################
  print('# Calculate C-domain/C-helix angle (ang_CHs)...')
  start = time.perf_counter()
  m_df.ang_CHs = VectorAngle([ d_df.hlx_cent.to_numpy(), d_df.dfg_d_ca.to_numpy(),
                               d_df.hlx_cent.to_numpy(), d_df.c_glu_cg.to_numpy() ])
  end = time.perf_counter()
  print(' {0:.1f} ms for {1} frames\n'.format((end-start)*1000, frames))
  
####################
  print('# Calculate N-/C-domain distance (dist_NC)...')
  start = time.perf_counter()
  m_df.dist_NC = Distance([ d_df.b3k_ca.to_numpy(), d_df.dfg_d_ca.to_numpy() ])
  end = time.perf_counter()
  print(' {0:.1f} ms for {1} frames\n'.format((end-start)*1000, frames))

####################
  print('# Calculate N-domain/C-helix distance (dist_NH)...')
  start = time.perf_counter()
  m_df.dist_NH = Distance([ d_df.b3k_cg.to_numpy(), d_df.c_glu_cg.to_numpy() ])
  end = time.perf_counter()
  print(' {0:.1f} ms for {1} frames\n'.format((end-start)*1000, frames))

####################
  print('# Calculate C-domain/C-helix distance (dist_CH)...')
  start = time.perf_counter()
  m_df.dist_CH = Distance([ d_df.dfg_d_ca.to_numpy(), d_df.c_glu_ca.to_numpy() ])
  end = time.perf_counter()
  print(' {0:.1f} ms for {1} frames\n'.format((end-start)*1000, frames))

  return m_df


##########################################################################
## read in trajectory using the supplied topology file (template pdb)
class ReadTraj(object):
  def __init__( self, top='' ):
    self.top = top

  def __call__( self, traj_file ):
    return md.load( traj_file, top=self.top )

#########
## each elements in this object will contain a list of coordinates
class TrajCoord(object):
  def __init__( self, ca=[], cb=[], cg=[], cn=[] ):
    self.ca = ca
    self.cb = cb
    self.cg = cg
    self.cn = cn

##########################################################################
## return an array of user-defined residue atoms, also check if residue is Gly,
## use average of HA2 and HA3 coords as CB and CG.
## selection with MDtraj, similar to VMD selection syntax
## For the collected coordinates, somehow MDTraj is missing by 1 decimal place
## and need to multiple the actual coords by a factor of 10
def SelectAtom( traj, resid, around=0, pkl=False ):

  Coord = TrajCoord()    # create object to hold atom coords

  ## if true, include extra residues before and after "resid"
  if around:
    resid = '{0} to {1}'.format(resid-around, resid+around)
  else:
    # md.Trajectory.xyz[frames, selection, atom-coordinates]
    select_ca = traj.top.select('resi {0} and name {1}'.format(resid, 'CA'))
    Coord.ca  = np.array(list(zip(*(10*traj.xyz[:, select_ca, :])))[0] )
#    print('ca', select_ca,'n',Coord.ca)

  ## check if residue is Glycine, single residue, or several residues
  if CheckResidue(traj.topology, resid) == 'helix':
    ## input is a set of 7-9 residues of C-helix centering on conserved Glu
    ## to calculate the helix axis curve, use center-most atoms as coordinates
    select_bb = traj.top.select('(resi {0}) and (name CA C N)'.format(resid))
    mid_atom  = ArrayCent(len(select_bb)) # center-most atom of helix atoms
 #   print('bb', select_bb)
    Frames    = 10*np.array(traj.xyz[:, select_bb, :], dtype=np.float64)
    frames = len(Frames)

#################
    ## 2nd-order regression of C-helix, take mid-point as helix center coords
    ## ** Using chunks of 'apply' + MPI, it is no different than MPI alone
    ##    but faster than 'loops' and 'apply' alone
    ## Use map will improve performance over imap by ~8-10% due to mpi overhead
    ## 2nd-order regression of C-helix, take mid-point as helix center coords
    start = time.perf_counter()
    print('# Calculate 2nd-order regression on C-helix for axis mid-point...')
    mpi       = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    Reg2_list = [x for x in tqdm(mpi.map(CalculateHelixAxis,Frames),total=frames)]
    end = time.perf_counter()
    print(' {0:.1f} ms for {1} frames with MPI\n'.format((end-start)*1000, frames))

#    def CalculateHelixAxis_InParallel( chunk ):
#      chunk_df = pd.DataFrame( {'frame': list(chunk)} )
#      Reg2_list = chunk_df.apply(lambda row: CalculateHelixAxis(row['frame']), axis=1 )
#      return Reg2_list

#    start = time.perf_counter()
#    f_chunks = np.array_split(Frames, multiprocessing.cpu_count())
#    print('# Calculate 2nd-order regression on C-helix for axis mid-point...')
#    mpi        = multiprocessing.Pool(processes=multiprocessing.cpu_count())
#    Chunk_list = [x for x in tqdm(mpi.imap(CalculateHelixAxis_InParallel, f_chunks),total=len(f_chunks))]
#    Reg2_list = [item for sublist in Chunk_list for item in sublist]
#    end = time.perf_counter()
#    print(' {0:.1f} ms for {1} frames with MPI+Apply\n'.format((end-start)*1000, frames))
 
#    start = time.perf_counter()
#    df['reg2_a'] = df.apply(lambda row: CalculateHelixAxis(row['frame']), axis=1 )
#    Reg2_list = df['reg2_a'].to_list()
#    end = time.perf_counter()
#    print(' {0:.1f} ms for {1} frames with Apply\n'.format((end-start)*1000, frames))

#################     
    Coord.cn  = np.asarray([Reg2[mid_atom] for Reg2 in Reg2_list])
    mpi.close()
    mpi.join()

  ## most residues have 'CB' and 'CG', but some have branched 'CG*' and 'CD*'
  ## Branched: Val/Ile have 'CG*', Cys has 'SG', Thr has 'OG*|CG*'
  ## Met has 'SD', Asp/Asn have 'OD*|ND*', His/Phe/Tyr/Trp have 'CD*|ND*'
  ## if CD not available, use CG; if CG not available, use CB
  elif CheckResidue(traj.topology, resid) != 'GLY':
    select_cb = traj.top.select('resi {0} and name {1}'.format(resid,'CB'))
#    print('cb', select_cb)
    Coord.cb  = np.array(list(zip(*(10*traj.xyz[:, select_cb, :])))[0], dtype=np.float64)

    select_cg = traj.top.select('resi {0} and (name =~ "{1}")'.format(
                                 resid, 'CG|OG|SG'))
    if len(select_cg) == 1:
      Coord.cg = np.array(list(zip(*(10*traj.xyz[:, select_cg, :])))[0], dtype=np.float64)
    elif len(select_cg) > 1:
      Frames   = 10*np.array(traj.xyz[:, select_cg, :], dtype=np.float64)
      Coord.cg = np.asarray(list(zip(*[np.mean(frame, axis=0) for frame in Frames]))[0])
    else:
      Coord.cg = Coord.cb
#    print('cg', select_cg, '\n', Coord.cg)

  else:   # if it is Glycine, use HA2 HA3 average as substitute of CB and CG
    topx, bonds = traj.topology.to_dataframe()
    print( ( topx[topx['resSeq']==resid+1] ) )
    select_h = traj.top.select('resi {0} and (name =~ "{1}")'.format(resid,'HA'))
    Frames   = 10*np.array(traj.xyz[:, select_h, :], dtype=np.float64)
    AvgTraj  = np.array([ list(np.mean(frame, axis=0)) for frame in Frames ])
    Coord.cb = AvgTraj
    Coord.cg = AvgTraj
    print(' * GLY{0} found, use HA* as CB/CG coords'.format(resid+1))
#    print('h', select_h, '\n', Coord.cg)
  return Coord


##########################################################################
## the numbers in topology record is same as user input (1-based), but actual
## mdtraj internal record (0-based) is 1 fewer than user input
def CheckResidue( top, resid ):

  ## if input is a set of residue, skip
  if re.search('to', str(resid)):
    print(resid, 'in GLU')
    return 'helix'

  top_df, bonds = top.to_dataframe()
#  print(top_df)
  print(list(top_df[top_df['resSeq'] == resid+1]['resName'])[0], resid+1)

  if list(top_df[top_df['resSeq'] == resid+1]['resName'])[0] == 'GLY':
    return 'GLY'
  else:
    return 'other'


##########################################################################
class CollectCoords(object):
  def __init__( self, coord_dfg_d=[], coord_dfg_f=[], coord_b3k=[],
                      coord_c_glu=[], coord_hlx=[] ):
    self.coord_dfg_d = coord_dfg_d
    self.coord_dfg_f = coord_dfg_f
    self.coord_b3k   = coord_b3k
    self.coord_c_glu = coord_c_glu
    self.coord_hlx   = coord_hlx

##########################################################################

class ExtractCoords(object):
  def __init__( self, dfg='', b3k='', c_glu='', pkl=False ):
    self.dfg   = dfg
    self.b3k   = b3k
    self.c_glu = c_glu
    self.pkl   = pkl

##################
  def __call__( self, inp ):
    return self._extract_coords( inp )
    
  def _extract_coords( self, inp ):
    Crd = CollectCoords()

    ## extract traj-coordinates of each kinase residues, mdtraj uses 0-based ID
    Crd.coord_dfg_d = SelectAtom(inp, int(self.dfg)-1,   around=0)
    Crd.coord_dfg_f = SelectAtom(inp, int(self.dfg)+0,   around=0)
    Crd.coord_b3k   = SelectAtom(inp, int(self.b3k)-1,   around=0)
    Crd.coord_c_glu = SelectAtom(inp, int(self.c_glu)-1, around=0)
    Crd.coord_hlx   = SelectAtom(inp, int(self.c_glu)-1, around=4)

    if self.pkl and not os.path.isfile(self.pkl):
      with open(self.pkl, 'wb') as fo:
        pickle.dump(Crd, fo, protocol=pickle.HIGHEST_PROTOCOL)
        print('  ## INFO: Write structural coords to: {0}\n'.format(self.pkl))

    return Crd

##########################################################################
## Calculate the helix axis using coordinates supplied, calculate 2nd-order
## regression curves to represent helix axis.
def CalculateHelixAxis( Coords ):
#  print('Coords\n', Coords)
  count = len(Coords)
#  print(count)

  # Linear regression on Cartesian Coordinates to estimate helix axis coords.
  # Use moving sets of points on both end to average out regression error
  # iterate range(3) to calculate x,y,z coordinates separately.
  # Use minimium of 7 residues (> 21 points),
  if count >= 15:
    posit  = 6                # 3 atoms (N,C,CA) = 1 residue
  else:
    posit  = 1
  xcount = count - posit    # reduced number of points to do LSQ

  Fn2Pts = []
  for m in range(0,posit):
#    Fn2 = [LsqFit_1d(range(xcount), Coords[m:m-posit, x] ,2) for x in range(3)]
    Fn2 = LsqFit_nd(range(xcount), Coords[m:m-posit], 2)
    Fn2Pts.append( [np.asarray([f(x) for f in Fn2]) for x in range(count) ])
  Reg2 = np.mean(Fn2Pts, axis=0)
  return Reg2


##########################################################################
## Take in coordinates, calculate vectors among the coords, generate
## Cross-Products of the pairs
## Asp-CG (r1), Asp-CA (r2), Phe-CA (r3), Phe-CG (r4)
def CalculateDFGVectors( inp ):
  r1, r2, r3, r4 = inp
  x = len(r1)

  # need to convert the list of vectors into dict before into Dataframe
  vec = {
  'r21': r1 - r2,  # (AspCG - AspCA)      D (r1)
  'r23': r3 - r2,  # (AspCA - PheCA)       \______ (r3)
  'r32': r2 - r3,  # (PheCA - AspCA)      (r2)    \
  'r34': r4 - r3,  # (PheCG - PheCA)               F (r4)

  'temp1':  np.zeros(x), 'temp2':  np.zeros(x) }
  t_df = pd.DataFrame(vec, index=range(x))

  ur21 = t_df.r21.div( VecMag( list(t_df.r21) ) )   # univector of r21
  ur23 = t_df.r23.div( VecMag( list(t_df.r23) ) )
  ur32 = t_df.r32.div( VecMag( list(t_df.r32) ) )
  ur34 = t_df.r34.div( VecMag( list(t_df.r34) ) )

  ## To enable pandas vectorized calculation, transpose the Mx3 to 3xM format
  ## for calculation, which is later transposed back to Mx3 format when done
  t1 = {'temp1': list(VecCross( [np.array(list(ur21)).T, np.array(list(ur23)).T] ).T) }
  t2 = {'temp2': list(VecCross( [np.array(list(ur34)).T, np.array(list(ur32)).T] ).T) }
  t_df.temp1 = pd.DataFrame(t1)
  t_df.temp2 = pd.DataFrame(t2)

  y = {'p1': t_df.temp1.div( VecMag(list(t_df.temp1)) ),  # univector for p1
       'p2': t_df.temp2.div( VecMag(list(t_df.temp2)) ),
       'v3': ur23 }                               # already univector for v3
  u_df = pd.DataFrame(y)
  del t_df

  # V = vector of the side chain, P = cross product vector of sidechain/D:D+1
#  p1 = np.cross( r21/VecMag(r21), r23/VecMag(r23) )
#  p2 = np.cross( r34/VecMag(r34), r32/VecMag(r32) )
#  v3 = r23/VecMag(r23)

  return u_df


##########################################################################
## Magnitude of a vector, equivalent to np.linalg.norm( v )
## this algorithm is much faster than np.linalg.norm function
def VecMag( v ):
  return np.sqrt((np.asarray(v)**2).sum(-1))

#################
## Distance between 2 points for pandas vectorization
def Distance( inp ):
  a, b = inp
  return VecMag( list( VecGen([a, b]) ) )

#################
# ## Generate a unit vector from 2 coordinates for pandas vectorization
def VecGen( inp ):
  a, b = inp
  return ( np.array(b) - np.array(a) )

#################
## Cross product in the most basic form for pandas vectorization
def VecCross( inp ):
  a, b = inp
  c = np.array([ a[1]*b[2] - a[2]*b[1],
                 a[2]*b[0] - a[0]*b[2],
                 a[0]*b[1] - a[1]*b[0]  ] )
  return c

#################
## Dot product in the most basic form for pandas vectorization
def VecDot( inp ):
  a, b = inp
#  print('a\n', a, '\n', 'b', '\n', b)
  c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
  return c

#################
## Angle between vectors for pandas vectorization
def VectorAngle( inp ):
  r1, r2, r3, r4 = inp
  
  o_df = pd.DataFrame(index=range(len(r1)), columns=['v1', 'v2', 'ang'])

  o_df.v1 = ( VecGen( [r1, r2] ) )   # df of v1 vector
  o_df.v2 = ( VecGen( [r3, r4] ) )
#  print('vecgen L592\n', r1, '\n', r2, '\n', o_df.v1)

  uv1 = o_df.v1.div( VecMag( list(o_df.v1)) )   # df of v1 univector
  uv2 = o_df.v2.div( VecMag( list(o_df.v2)) )
  
  ## To enable pandas vectorized calculation, transpose the Mx3 to 3xM format
  ## for calculation, which is later transposed back to Mx3 format when done
  dot = VecDot([ np.array(list(uv1)).T, np.array(list(uv2)).T ]).transpose()
  o_df.ang = np.arccos( dot ) * 180/np.pi
#  print('o_df L602\n', o_df)
  return o_df.ang

#################
## n-order Least-square fit of 2 arrays, x and y, return an object/def as
## a nst-order polynomial function f(x) = mx+c, or f(x) = mx^2+nx+c, etc
## i.e.                            f(1) = 15       f(5) = 25
## polyfit of single nd-array seems slightly faster than polyfit of multiple
## 1d-array, probably due to numpy's vectorization
def LsqFit_1d( X, Y, order=1 ):
  return np.poly1d( np.polyfit( X, Y, order, full=False)) # for single polyfit
def LsqFit_nd( X, Y, order=1 ):
  Fits = list( zip(*np.polyfit( X, Y, order, full=False )) )
  return [ np.poly1d(coeff) for coeff in Fits ]           # for multiple polyfit

#################
# Find the center element in an array. The result is 1 less than correct
# answer -- Python numbering starts from 0. i.e. 9-element array is [0,...,8]
# center is 4, not 5
def ArrayCent( count ):
  if count % 2 == 0:
    center = (count/2)-1
  else:
    center = ((count-1)/2)
  return int(center)


##########################################################################
def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-templ', dest='tmpl_file', required=False,
                 help='Template PDB structure (exact match to Topology Atom List and aligned to Ref structure 1ATP)')
  p.add_argument('-traj', dest='traj_file', required=False,
                 help='Trajectory file, or an ordered list of traj filenames')
  p.add_argument('-out', dest='outpref', required=False,
                 help='Output prefix')

  p.add_argument('-pkl', dest='pkl', required=False,
                 help='Use pre-pickled data (def: False)')
  p.add_argument('-b3k', dest='b3k', required=False,
                 help='Residue Number of (beta-3 Lys)')
  p.add_argument('-dfg', dest='dfg', required=False,
                 help='Residue Number of (DFG Asp)')
  p.add_argument('-glu', dest= 'c_glu', required=False,
                 help='Residue Number of (C-helix Glu)')

  p.add_argument('-superp', dest= 'superp', required=False,
                 help='*Optional: VMD-like selection string to perform superposition (default: False)')

  args=p.parse_args()
  return args

##########################################################################
if __name__ == '__main__':
  main()

##########################################################################
#
# Peter M.U. Ung @ MSSM/Yale
#
# v1.0  19.05.06  finished reading and calculating traj coodinates
#                 optimized with DataFrame/Numpy vectorization and MPI
# v2.0  19.05.18  implement R RandomForest classifier to assign conformation
#
#
