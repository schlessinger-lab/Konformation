#!/usr/bin/python

##########################################################################
#
# v 1.0 - 17.08.17  # modified from 0_check_type_I-II.py
#   
###
### ignore this script, turns out it is useless, classifying things wrong :p
###
#
# distance between the center-of-mass (COM) of the inhibitor to the COM 
# of a known type-I (ligand 215 from 2FB8 for S/T-kinases and XIN from 3C7Q
# for Y-kinases) or type-II ( imitanib STI from 3HEC and 2PL0 for S/T- and 
# Y-kinases, respectively)  kinase inhibitor was calculated. 
# Ligands with COMtype-I > 5.00 A and COMtype-II < 3.05 A are considered 
# type-II ligands (118 ligands), while all others are considered type-I ligands
# (1340 ligands). In addition to these ligands, analogs of the type-II inhibitor
# PD5 from 3EL81 (AD-36, AD-57, AD-58, AD-80, and AD-81) that have wide 
# anti-kinase spectrum was added to the ligand set.
#
#   S/T/D-kinases:
#     [ 2FB8_A.1atp.pdb 215 ] + [ 3HEC_A.1atp.pdb STI ]
#
#   Y-kinases:
#     [ 3C7Q_A.2bdf.pdb XIN ] + [ 2PL0_A.2bdf.pdb STI ]
#
#   MUST REVIEW the ligands individually, especially Type_I-B and those with
#   COM(type-II) and COM(type-I) very close in value
#
##########################################################################

import sys
msg = '''
    > {0}
        [ List of PDB ] [ PDB Directory ]
        [ Type-I  Sample PDB ] 
        [ Type-II Sample PDB ]
        [ Cutoff Distance ]
        [ Output file ]\n
      * STD_kinase: 2FB8_A.1atp.pdb 215 3HEC_A.1atp.pdb STI | 3.30 
        Y_kinase:   3C7Q_A.2bdf.pdb XIN 2PL0_A.2bdf.pdb STI | 3.50\n
#   MUST REVIEW the ligands individually, especially Type_I-B and those with
#   COM(type-II) and COM(type-I) very close in value\n
#   Also watch out the ones many atoms\n'''.format(sys.argv[0])
if len(sys.argv) != 7: sys.exit(msg)

import os,re,glob
import numpy as np
from CommonUtility import *
import multiprocessing
from aa_residue import *

##########################################################################
def main(pdb_list, directory, typeI_pdb, 
         typeII_pdb, cutoff, outfile):

  # read in the PDB list
  with open(pdb_list, 'r') as fi:
    PDB = [l.rstrip() for l in fi]

  # Read thru all PDBs in the directory to get the PDB_ID of those Ligands
  PDB_Ligs   = CollectHetero(directory, PDB)
  TypeI_Lig  = CollectHetero(directory, [typeI_pdb])
  TypeII_Lig = CollectHetero(directory, [typeII_pdb])

  TypeI_COM  = CalculateCOM(TypeI_Lig[0] )
  TypeII_COM = CalculateCOM(TypeII_Lig[0])

  Dists = []
  Data  = []
  for Itm in PDB_Ligs:
    pdb = Itm[0]
    if len(Itm[1]) == 0:
      Data.append( [pdb, '-', '-', 0, 0., 0., 'NaN'] )
      continue

    Lig_COM = CalculateCOM(Itm)
    lig     = Lig_COM[0][1]
    ha_num  = Lig_COM[0][3]
    resi    = Lig_COM[0][4]

    dist_I  = np.linalg.norm(TypeI_COM[0][2]  - Lig_COM[0][2])
    dist_II = np.linalg.norm(TypeII_COM[0][2] - Lig_COM[0][2])
    
#    if dist_I + dist_II < 22.0:  # ligand within binding site
#      if   dist_II <= cutoff:
#        lig_type = 'Type_II'    # ok
#      else:
#        lig_type = 'Type_I'
#    else:
#      lig_type = 'other'        # ok

    if dist_I + dist_II < 22.0:  # ligand within binding site
      if dist_I > cutoff*2:
        if dist_II < cutoff:
          lig_type = 'Type_II'
        elif dist_II > cutoff*2:
          lig_type = 'Type_I'
        else:
          lig_type = 'Type_II'

      elif dist_I > cutoff:
        if dist_II < cutoff:
          lig_type = 'Type_II'
        elif dist_II > cutoff*2:
          lig_type = 'Type_I'
        else:
          lig_type = 'Type_I.5'

      else:     # dist_I <= cutoff
        if dist_II > cutoff*2:
          lig_type = 'Type_I'
        elif dist_II > cutoff:
          lig_type = 'Type_I.5'
        else:   # dist_II <= cutoff
          lig_type = 'Type_II'
    else:
      lig_type = 'Other'

    Data.append( [pdb, lig, resi, ha_num, dist_I, dist_II, lig_type] )


  Rslt = sorted( Data, key=lambda x: (x[6], x[3]) )
  with open(outfile, 'w') as fo:
    fo.write('PDB\tlig\tresi\tATM\t \tdist-I\tdist-II\t \tlig-type\n')
    for l in Rslt:
      fo.write('{0}\t{1:3s}\t{2:4s}\t{3:3d}\t=\t{4:6.3f}\t{5:6.3f}\t=\t{6}\n'.format(
                l[0], l[1], l[2], l[3], l[4], l[5], l[6]) )


################################
class ExtractHetero(object):
  def __init__(self, directory=''):
    self.dir = directory

  def __call__(self, pdb_name):
    return self.extract_pdb(pdb_name)

#######################
  def extract_pdb(self, pdb_name):
    try:
      os.path.isfile(self.dir+'/'+pdb_name)
    except IOError:
      print('  Error: Cannot find file: '+self.dir+'/'+pdb_name)
      return None

    with open('{0}/{1}'.format(self.dir, pdb_name), 'r') as fi:
      Lines = [ l for l in fi if re.search(r'^HETATM', l) ]
    pdb_id = pdb_name.split('.')[0]
    return [pdb_id, Lines]

################################

##########################################################################
def CollectHetero(directory, PDB):

  mpi = multiprocessing.Pool(processes = multiprocessing.cpu_count())
  Het = ExtractHetero(directory=directory)
#  Tmp = [Het(pdb) for pdb in PDB]
  Tmp = mpi.map(Het, PDB)
  mpi.close()
  mpi.join()

  Data = []
  for Itm in Tmp:
    if Itm is not None:
      Hetero = {}
      if len(Itm[1]) > 0:
        for l in Itm[1]:
          lig = l[17:20]            # molecule res_nam
          mol = l[22:26]            # molecule res_id
          if SaltAdditive(lig):     # remove additives, salts, HOH
            pass
          else:
            if lig in Hetero:
              if mol in Hetero[lig]:
                Hetero[lig][mol].append(l)
              else:
                Hetero[lig][mol] = [l]
            else:
              Resi = { mol: [l] }
              Hetero[lig] = Resi
    Data.append( [Itm[0], Hetero] )

  return Data


##########################################################################
def CalculateCOM( Input ):
  pdb_id, Hetero = Input
  
  Rslt = []
  if len(Hetero) > 0:
    for lig in Hetero:
      Coord = [[0,0,0],0]     # [Coordinates, no. of atom]

      for mol in Hetero[lig]:
        for line in Hetero[lig][mol]:
          Coord[0][0] += float(line[30:38])
          Coord[0][1] += float(line[38:46])
          Coord[0][2] += float(line[46:54])
          Coord[1]    += 1
        Rslt.append( [pdb_id, lig, np.divide(Coord[0], Coord[1]),
                      Coord[1], mol] )
  else:
    Rslt.append( [pdb_id, None, None, None, None] )
    
  return Rslt
  

##########################################################################
if __name__ == '__main__':
  main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], 
       float(sys.argv[5]), sys.argv[6])

