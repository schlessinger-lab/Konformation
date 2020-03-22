#!/usr/bin/env python3

import sys
##########################################################################
#
#   Peter M.U. Ung  @ MSSM
#
#   ** KinMetrics -- calculate descriptors for kinase conformation 
#                    classification (Kinformation)
#                    (Optional) check input PDB by getting its sequence,
#                    align to kinome to remove non-kinases
#                    (Optional) align input PDB to reference 1ATP base
#
#   Purpose:  parse input PDB, make sure they are aligned to 1ATP base. 
#       Running profile alignment to align PDB sequence to 1ATP sequence,
#       get the residue positions in sequence to calculate descriptors,
#       supply descriptors to R script RandomForest to determine conformation
# 
##########################################################################

from x_konf_vars import GenerateTemplSetupScript

msg = '''
  > {0}
      <file>   [ Parameter File ]\n
  Optional:
      -t       [ Generate template input file ]
'''.format(sys.argv[0])
if len(sys.argv) == 1: sys.exit(msg)
if sys.argv[1] == '-t':
  templ_file = 'Template_parameter_file.in'
  GenerateTemplSetupScript(templ_file)  # Generate template param file
  sys.exit('\n** Generated - \033[31m{0}\033[0m\n'.format(templ_file))

import os,re
import time
import numpy as np
import pandas as pd

from Bio import SeqIO
from tqdm import tqdm
from shutil import copy2
from pathos import multiprocessing
from sklearn.impute import SimpleImputer

from x_konformation  import Konformation
from x_pdb_superpose import SuperposePDB
from x_search_align  import CacheSeqDatabase
from x_konf_vars     import ParseParameterFile
from x_search_align  import RemoveFastaGapColumn
from x_search_align  import CheckInputStructures
from x_search_align  import GenerateProfileAlignment

from x_kinfo_SK_classify import KinfoClassify

##########################################################################

def main( param_list ):

  # Read in parameters
  parm   = ParseParameterFile(param_list)
  outext = parm['OUTEXT'][0]
  output_prefix = parm['OUTPREF'][0]

  hom_dir = parm['HOMEDIR'][0]
  pdb_dir = parm['PDBDIR'][0]
  rst_dir = hom_dir+'/1_result/'
  tmp_dir = hom_dir+'/x_work/'

  print('  \033[34m> Script directory:\033[0m '+parm['SCRIPT'][0])
  print('  \033[34m> Internal PDB directory:\033[0m '+parm['PDBDIR'][0])
  print('  \033[34m> Reference kinase PDB:\033[0m '+parm['REFPDB'][0])
  print('  \033[34m> Reference kinase residues:\033[0m '+parm['REFRES'][0])
  print('  \033[34m> Output data prefix:\033[0m '+parm['OUTPREF'][0])

  if not os.path.exists(rst_dir):
    os.makedirs(rst_dir)
  if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)

  os.chdir(tmp_dir)       # Go to temporary directory to start working
  
  refname = parm['REFPDB'][0].split('/')[-1].split('.pdb')[0]
  Fasta_Dict = CacheSeqDatabase(parm['FASTA'][0])


#####################################
  # read Target PDB from list, check if structure has multiple chains
  # regenerate PDB and repopulate target PDB list with chain_id
  if parm['CHECKPDB'][0]:
    print('\033[31m  ## CHECKPDB: Input PDB List have NOT been checked for kinase domain. Checking... ##\033[0m')
    Temp1 = pd.read_csv(parm['PDBLIST'][0], sep='\s+', header=None, comment='#').iloc[:,0].values.tolist()

    if parm['MPICPU'][0] == 1:
      Temp2 = [CheckInputStructures(pdb) for pdb in tqdm(Temp1)]
    else:
      if parm['MPICPU'][0] == 0:
        mpi_cpu = multiprocessing.cpu_count()
      else:
        mpi_cpu = parm['MPICPU'][0]
      mpi   = multiprocessing.Pool(mpi_cpu)
      Temp2 = [x for x in tqdm(mpi.imap(CheckInputStructures,Temp1),total=len(Temp1))]
      mpi.close()
      mpi.join()

    Targets = [item for sublist in Temp2 for item in sublist]
  else:
    print('\033[34m  ## CHECKPDB: Input PDB have been Pre-Checked for kinase domain ##\033[0m')
    print(os.getcwd())
    Temp = pd.read_csv(parm['PDBLIST'][0], header=None, comment='#', sep='\s+').iloc[:,0].to_numpy()
    Targets = [ pdb_dir+'/'+pdb for pdb in Temp ]

  print('\n ** Input Target PDB: \033[31m{0}\033[0m'.format(len(Targets)))


#####################################
  # check if PDB has been superposed to 1ATP; if not, superpose and rename
  if parm['SUPERPOSED'][0]:
    print('\n\033[34m  ## SUPERPOSED: Input PDB have been Pre-Superposed to reference kinase ##\033[0m')
    for pdb in tqdm(Targets):
      if not os.path.isfile(pdb):
        copy2(pdb, rst_dir)
  else:
    print('\n\033[31m  ## SUPERPOSED: Input PDB have NOT been Pre-Superposed to reference kinase. Superposing... ##\033[0m')
    SuperposePDB( parm['PYMOL'][0], parm['REFPDB'][0], Targets, 'pdb', 
                  outext, parm['REFRES'][0], output_prefix )
    Temps = []
    for pdb in tqdm(Targets):
      Temps.append('{0}/{1}.{2}'.format(tmp_dir,pdb.split('/')[-1].split('.')[0],outext))
      copy2('{0}/{1}.{2}'.format(tmp_dir,pdb.split('/')[-1].split('.')[0],outext), rst_dir)
#      Temps.append(tmp_dir+'/'+pdb.split('/')[-1].split('.')[0]+'.'+outext)
#      copy2(tmp_dir+'/'+pdb.split('/')[-1].split('.')[0]+'.'+outext, rst_dir)
    Targets = Temps


#####################################
  # Get FASTA from PDB structures, find most similar seq, do Profile align,
  # reject input sequence with identity match < 30% to known kinase structures
  ## For new PDB, need to extract its FASTA, use blastdb to find 10 closest FASTA MSA,
  ## align to MSA to create correct gapping, then use the newly aligned FASTA to find
  ## the correct key residues for extraction
  if parm['CHKALIGN'][0]:
    print('\n\033[31m  ## CHKALIGN: Need to check PDB for non-kinases and generate MSA Alignment. Checking... ##\033[0m')
    oAlign = GenerateProfileAlignment( 
              hom_dir = hom_dir,  rst_dir = rst_dir,  
              tmp_dir = tmp_dir,  ref_pdb = refname,
              f_nogap = parm['BLASTDB'][0], f_dict  = Fasta_Dict )

    if parm['MPICPU'][0] == 1:
      Temps = [ oAlign(pdb) for pdb in tqdm(Targets) ]
    else:
      if parm['MPICPU'][0] == 0:
        mpi = multiprocessing.Pool()
      else:
        mpi   = multiprocessing.Pool(parm['MPICPU'][0])
      Temps = [x for x in tqdm(mpi.imap(oAlign, Targets),total=len(Targets))]
      mpi.close()
      mpi.join()

  # Remove any input kinase that have been rejected due to low seq identity
    with open('_TEMP.temp_comb.fasta', 'w') as fo:
      SeqIO.write( [Fasta_Dict[refname]]+[x for x in Temps if x is not None],
                    fo, 'fasta' )

    # Remove gaps in temp fasta file
    RemoveFastaGapColumn('_TEMP.temp_comb.fasta', '{0}.all_comb.fasta'.format(output_prefix))

    # Reset MSA-aligned Fasta file in working directory for use in Kinformation
    parm['FASTA'][0] = tmp_dir+'{0}.all_comb.fasta'.format(output_prefix)

    # If input PDB have to be superposed, new input PDB will be in tmp_dir, reset
    parm['PDBDIR'][0] = tmp_dir
  else:
    print('\n\033[34m  ## CHKALIGN: PDB have been Pre-Checked for non-kinases and Generated MSA Alignment ##\033[0m')
    parm['FASTA'][0] = parm['PDBALIGN'][0]

  ## Move to result directory to use data there to calculate conformation
  os.chdir(rst_dir)


###################

  ## Run Konformation to examine kinase conformation
  print('\n\033[34m  ## Run Konformation to calculate structural parameters ##\033[0m')
  parm_df = Konformation( parm, Targets, output_prefix )

  print('\n#########################################\n')

##################
  ## Columns used for Conformation Classification
  Selected_Cols = [ 'p1p1x','p2p2x','r3r3x','h_cgvc',
                    'ang_NHs','ang_CHs','dist_NH','dist_CH' ]
  df = parm_df[Selected_Cols]

  ## Impute any missing value in dataset before ML classication
  imputer = SimpleImputer(missing_values=np.nan, strategy='median')
  Imp_lst = imputer.fit_transform( df.to_numpy() )  # skip 1st column, pdb_id

  data_df = pd.DataFrame(Imp_lst, index=df.index, columns=df.columns)
  data_df.reset_index(inplace=True)

#################
  ## use Kinformation Random Forest Classifier to assign conformation/confidence
  KinfoClassify(  data_df, parm['DATADIR'][0], output_prefix, parm,
                  parm['USESKL'][0]  )

  print('\n#########################################\n')


#########################################################################
if __name__ == '__main__':
  main( sys.argv[1] )

##########################################################################
#
#   v1.0  18.03.11
#   v2.0  20.01.12
#   v3.0  20.02.27  add KinfoClassify function to do conformation classification
#
