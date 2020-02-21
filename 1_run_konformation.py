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

templ_file = 'Template_parameter_file.in'
msg = '''\n\t> {0}\n\t\t[Parameter File]
\t** Generated - \033[31m{1}\033[0m\n'''.format(sys.argv[0], templ_file)
if len(sys.argv) != 2: 
  GenerateTemplSetupScript(templ_file)  # Generate template param file
  sys.exit(msg)

import os,re
import pandas as pd

from Bio import SeqIO
from tqdm import tqdm
from shutil import copy2
from pathos import multiprocessing

from x_konformation  import Konformation
from x_pdb_superpose import SuperposePDB
from x_search_align  import CacheSeqDatabase
from x_konf_vars     import ParseParameterFile
from x_search_align  import RemoveFastaGapColumn
from x_search_align  import CheckInputStructures
from x_search_align  import GenerateProfileAlignment

##########################################################################

def KinDentity( param_list ):
  
  # Read in parameters
  RefRes = ParseParameterFile(param_list)
  outext = RefRes['OUTEXT']
  output_prefix = RefRes['OUTPREF']

  hom_dir = RefRes['HOMEDIR']
  rst_dir = hom_dir+'/1_result'
  tmp_dir = hom_dir+'/x_work'

  print('  \033[34m> Script directory:\033[0m '+RefRes['SCRIPT'])
  print('  \033[34m> Internal PDB directory:\033[0m '+RefRes['PDBDIR'])
  print('  \033[34m> Reference kinase PDB:\033[0m '+RefRes['REFPDB'])
  print('  \033[34m> Reference kinase residues:\033[0m '+RefRes['REFRES'])
  print('  \033[34m> Output data prefix:\033[0m '+RefRes['OUTPREF'])

  if not os.path.exists(rst_dir):
    os.makedirs(rst_dir)
  if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)

  os.chdir(tmp_dir)       # Go to temporary directory to start working
  
  refname = RefRes['REFPDB'].split('/')[-1].split('.pdb')[0]
  Fasta_Dict = CacheSeqDatabase(RefRes['FASTA'])


#####################################
  # read Target PDB from list, check if structure has multiple chains
  # regenerate PDB and repopulate target PDB list with chain_id
  if RefRes['CHKPDB']:
    print('\033[31m  ## CHKPDB: Input PDB have NOT been checked for kinase domain. Checking... ##\033[0m')
    Temp1 = pd.read_csv(hom_dir+'/'+RefRes['PDBLIST'], sep='\s+', header=None, comment='#').iloc[:,0].values.tolist()

    mpi   = multiprocessing.Pool()
    Temp2 = [x for x in tqdm(mpi.imap(CheckInputStructures,Temp1),total=len(Temp1))]
#    Temp2 = [CheckInputStructures(pdb) for pdb in tqdm(Temp1)]
    mpi.close()
    mpi.join()

    Targets = [item for sublist in Temp2 for item in sublist]
  else:
    print('\033[34m  ## CHKPDB: Input PDB have been Pre-Checked for kinase domain ##\033[0m')
    Temp = pd.read_csv(hom_dir+'/'+RefRes['PDBLIST'], header=None, comment='#', sep='\s+').iloc[:,0].values.tolist()
    Targets = [ hom_dir+'/'+pdb for pdb in Temp ]

  print('\n ** Input Target PDB: \033[31m{}\033[0m'.format(len(Targets)))
  for pdb in Targets: print(pdb)


#####################################
  # check if PDB has been superposed to 1ATP; if not, superpose and rename
  if RefRes['SUPERPOSED']:
    print('\n\033[34m  ## SUPERPOSED: Input PDB have been Pre-Superposed to reference kinase ##\033[0m')
    for pdb in Targets:
      copy2(pdb, rst_dir)
  else:
    print('\n\033[31m  ## SUPERPOSED: Input PDB have NOT been Pre-Superposed to reference kinase. Superposing... ##\033[0m')
    SuperposePDB( RefRes['PYMOL'], RefRes['REFPDB'], Targets, 'pdb', 
                  outext, RefRes['REFRES'], output_prefix )
    Temps = []
    for pdb in Targets:
      Temps.append(tmp_dir+'/'+pdb.split('/')[-1].split('.')[0]+outext)
      copy2(tmp_dir+'/'+pdb.split('/')[-1].split('.')[0]+outext, rst_dir)
    Targets = Temps


#####################################
  # Get FASTA from PDB structures, find most similar seq, do Profile align,
  # reject input sequence with identity match < 30% to known kinase structures
  ## For new PDB, need to extract its FASTA, use blastdb to find 10 closest FASTA MSA,
  ## align to MSA to create correct gapping, then use the newly aligned FASTA to find
  ## the correct key residues for extraction

  if RefRes['CHKALIGN'] is None:
    print('\n\033[31m  ## CHKALIGN: Need to check PDB for non-kinases and generate MSA Alignment. Checking... ##\033[0m')
    oAlign = GenerateProfileAlignment( 
              hom_dir = hom_dir,  rst_dir = rst_dir,  
              tmp_dir = tmp_dir,  ref_pdb = refname,
              f_nogap = RefRes['BLASTDB'], f_dict  = Fasta_Dict )
    mpi = multiprocessing.Pool()
    Temps = [x for x in tqdm(mpi.imap(oAlign, Targets),total=len(Targets))]
#    Temps = [ oAlign(pdb) for pdb in tqdm(Targets) ]
    mpi.close()
    mpi.join()

  # Remove any input kinase that have been rejected due to low seq identity
    with open('_TEMP.temp_comb.fasta', 'w') as fo:
      SeqIO.write( [Fasta_Dict[refname]]+[x for x in Temps if x is not None],
                    fo, 'fasta' )

    # Remove gaps in temp fasta file
    RemoveFastaGapColumn('_TEMP.temp_comb.fasta', '{}.all_comb.fasta'.format(output_prefix))

    # Reset MSA-aligned Fasta file in working directory for use in Kinformation
    RefRes['FASTA'] = tmp_dir+'{}.all_comb.fasta'.format(output_prefix)

    # If input PDB have to be superposed, new input PDB will be in tmp_dir, reset
    RefRes['PDBDIR'] = tmp_dir
  else:
    print('\n\033[34m  ## CHKALIGN: PDB have been Pre-Checked for non-kinases and Generated MSA Alignment ##\033[0m')
    
  ## Move to result directory to use data there to calculate conformation
  os.chdir(rst_dir)



###################

  ## Run Konformation to examine kinase conformation
  print('\n\033[34m  ## Run Konformation to examine kinase conformation ##\033[0m')
  Konformation( RefRes, Targets, output_prefix )


#########################################################################
if __name__ == '__main__':
  KinDentity( sys.argv[1], sys.argv[2] )

##########################################################################
#
#   v1.0  18.03.11
#   v2.0  20.01.12
#
#
