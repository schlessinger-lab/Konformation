#!/usr/bin/python

import sys
##########################################################################
#
#   Peter M.U. Ung  @ MSSM
#
#   ** KinMetrics -- Process input PDB structures and calculate descriptors
#                    for kinase conformatoin classification (Kinformation)
#
#   Purpose:  parse input PDB, make sure they are aligned to 1ATP base. 
#       Running profile alignment to align PDB sequence to 1ATP sequence,
#       get the residue positions in sequence to calculate descriptors,
#       supply the descriptors to R script to determine kinase conformation
# 
##########################################################################

msg = '''\n\t> {0}\n\t\t[Parameter File]
\t\t[Output Prefix]\n'''.format(sys.argv[0])
if len(sys.argv) != 3: sys.exit(msg)

import os,re
from pathos import multiprocessing
from tqdm import tqdm
from shutil import copy2
from x_konformation  import *
from x_search_align  import *
from x_pdb_superpose import *
from x_konf_vars     import *

##########################################################################

def KinDentity( param_list, output_prefix ):
  
  # Read in parameters
  RefRes = ParseParameterFile(param_list)
  outext = RefRes['OUTEXT']

  hom_dir = RefRes['HOMEDIR']
  rst_dir = hom_dir+'/1_result'
  tmp_dir = hom_dir+'/x_work'
  if not os.path.exists(rst_dir):
    os.makedirs(rst_dir)
  if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)

  os.chdir(tmp_dir)       # Go to temporary directory
  
  refname = RefRes['REFPDB'].split('/')[-1].split('.pdb')[0]
  Fasta_Dic = CacheSeqDatabase(RefRes['FASTA'])


#####################################
  # read Target PDB from list, check if structure has multiple chains
  # regenerate PDB and repopulate target PDB list with chain_id
  with open(hom_dir+'/'+RefRes['PDBLIST'], 'r') as fi:
    mpi = multiprocessing.Pool(processes = multiprocessing.cpu_count())

    Temp1 = [hom_dir+'/'+pdb_name.rstrip() for pdb_name in fi]
#    Temp2 = [CheckInputStructures(pdb) for pdb in Temp1]
    Temp2 = [x for x in tqdm(mpi.imap_unordered(CheckInputStructures,Temp1),total=len(Temp1))]
    mpi.close()
    mpi.join()

    Targets = [item for sublist in Temp2 for item in sublist]
  print('\n ** Input Target PDB **')
  for pdb in Targets: print(pdb)


#####################################
  # check if PDB has been superposed to 1ATP; if not, superpose and rename
  if RefRes['SUPERPO']:
    for pdb in Targets:
      copy2(pdb, rst_dir)
  else:
    SuperposePDB( RefRes['REFPDB'], Targets, 'pdb', outext,
                  RefRes['REFRES'], RefRes['REMOVE'] )
    Temps = []
    for pdb in Targets:
      Temps.append(tmp_dir+'/'+pdb.split('/')[-1].split('.')[0]+outext)
      copy2(tmp_dir+'/'+pdb.split('/')[-1].split('.')[0]+outext, rst_dir)
    Targets = Temps


#####################################
  # Get FASTA from PDB structures, find most similar seq, do Profile align,
  # reject input sequence with identity match < 30% to known kinase structures
  mpi = multiprocessing.Pool(processes = multiprocessing.cpu_count())
  oAlign = GenerateProfileAlignment( 
           hom_dir = hom_dir,  f_bdb   = RefRes['BLASTDB'],
           ref_pdb = refname,  f_dic   = Fasta_Dic, 
           rst_dir = rst_dir,  tmp_dir = tmp_dir )
#  Temps = [ oAlign(pdb) for pdb in Targets ]
  Temps = [x for x in tqdm(mpi.imap_unordered(oAlign, Targets),total=len(Targets))]
  mpi.close()
  mpi.join()

  # Remove any input kinase that have been rejected due to low seq identity
  with open('_TEMP.temp_comb.fasta', 'w') as fo:
    SeqIO.write( [Fasta_Dic[refname]]+[x for x in Temps if x is not None],
                 fo, 'fasta' )

  os.system('t_coffee -other_pg seq_reformat -action +rm_gap 100 -in {0} > {1}'.format( '_TEMP.temp_comb.fasta', '_TEMP.all_comb.fasta' ))

  # Replace Fasta file in working directory for use in Kinformation
  RefRes['FASTA'] = tmp_dir+'/_TEMP.all_comb.fasta'

  # If input PDB have to be superposed, new input PDB will be in tmp_dir
  RefRes['PDBDIR'] = tmp_dir

  os.chdir(rst_dir)
  Konformation( RefRes, Targets, output_prefix )


#########################################################################
if __name__ == '__main__':
  KinDentity( sys.argv[1], sys.argv[2] )

##########################################################################
#
#   v1.0  18.03.11
#
#
#
