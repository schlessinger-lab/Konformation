#!/usr/bin/env python3

import sys,os,re
import pandas as pd

def DefaultVariables():

  pdb_dir = '/home/pmung/Dropbox/1_kinase/1_family/1_stdkinases/170109/2_align/'
  script  = '/home/pmung/Dropbox/9_scripts/3_program/structures/4_Konformation/'
  data_dir= script+'z_database/'

  parm = {
    'PDBDIR':   [pdb_dir,         '# Path to directory with PDB structures'],
    'SCRIPT':   [script,          '# Path to main scripts'],
    'DATADIR':  [data_dir,        '# Path to main scripts database'],
    'HOMEDIR':  [os.getcwd()+'/', '# Path to home directory for working'],
    'OUTPREF':  ['output_prefix', '# Output prefix for data generated'],

    'PDBLIST':  ['pdb.list', '# List of PDB to be examined; same as "WTNEWKNOWN" for 1_find_new_kinase_pdb.py'],
    'PDBALIGN': ['pdb.fasta','# Aligned FASTA of PDB corresponding to PDBLIST'],

    'CHECKPDB': ['False',    '# Initiate checking of input PDB provided in PDBLIST?'],
    'CHKALIGN': ['False',    '# Initiate checking of input PDB sequence to remove non-kinases?'],
    'SUPERPOSED':['True',   '# Has input PDB been pre-superposed to 1ATP ref residue?'],
    'MISSRES':  ['None',     '# User-input list of missing residue PDB, same directory to PDB'],

    'USERRF':   ['False', '# Use R::randomForest instead of SKLearn RFClassifier (def: False)'],
    'USESKL':   ['nn',    '# Use SKLearn ML model: rf|svm|nn|kn|dt|gp|gb (def: rf)'],

    'FASTA':    [data_dir+'MD_human_kinome_alignment.2019-2.fasta', '# FASTA file of aligned canonical human kinome, with gap'],
    'BLASTDB':  [data_dir+'MD_human_kinome_alignment.2019.nogap.fasta', '# FASTA file of unaligned canonical human kinome, NO gap'],
    'REFPDB':   [data_dir+'1ATP_E.pdb',   '# Reference PDB, bovine PKA (1ATP_E)'],
    'REFRES':   ['resi 121-139+162-183',  '# Reference 1ATP_E residues for PyMOL superposition'],
    'OUTEXT':   ['1atp.pdb',        '# Output extension for all 1atp-superposed PDB'],
    'PYMOL':    ['/usr/bin/pymol',  '# Full path to PyMOL Executable'],
    'MPICPU':   ['0',               '# Number of CPU to use (def: 0=all, 1-->any)'],


    'R_dfg_model':      'R_rf_model_dfg.190527.rda',
    'R_chx_model':      'R_rf_model_full.190527.rda',

    'kinfo_norm_param': 'kinfo_data_normalize_factor.171009.csv',

## 5-resid based reference sequences for searching other kinase PDB's relevant
## residue positions. Intended residue is always at the center of the sequence
## flanked by equal no. of residue on both sides. Need at least 5-resid to have
## good confidence in correct identification. Search algorithm will fail if PDB
## structure has missing or unnatural residue in corresponding sequence. Will
## need to manually find and save the missing residue into a MISSING PDB, and 
## append the info into the MISSING list. 
## For failed cases, try alternative 5-resid sequence that takes last residue

    'HELIX':    'HTLNEKRIL',  # [Ref 1ATP C-helix (Glu) | 5,7,9-resid center on E91]
    'NDOM':     'AMKIL',  # [Ref 1ATP N-lobe beta3(Lys)] K72
    'CDOM':     'VTDFG',  # [Ref 1ATP C-lobe DFG(Asp)] D184
    'DFGF':     'TDFGF',  # [Ref 1ATP C-lobe DFG(Phe)] F185

    'XHELIX':   'EKRIL',    # Alternative for HELIX, the first resid
    'ZNDOM':    'HYAMK',    # Alternative for NDOM, the last resid
    'ZCDOM':    'IQVTD',    # Alternative for CDOM, the last resid
    'ZDFGF':    'QVTDF',    # Alternative for DFG-F, the last resid

    'GATE':     'MVMEY',  # Gatekeeper M120
#    'RSPINE1':  'VKLEF',  # L106   R-spine of bovine PKA (1ATP_E)
#    'RSPINE2':  'RILQA',  # L95
#    'RSPINE3':  'TDFGF',  # F185
#    'RSPINE4':  'LIYRD',  # Y164
#    'CSPINE1':  'HYAMK',  # A70    C-spine of bovine PKA (1ATP_E)
#    'CSPINE2':  'NLLID',  # L173
#    'CSPINE3':  'ENLLI',  # L172
#    'CPSINE4':  'GEMFS',  # M128
#    'CSPINE5':  'YEMAA',  # M231
#    'CSPINE6':  'GRVML',  # V57
#    'CSPINE7':  'LLIDQ',  # I174
#    'CSPINE8':  'GVLIY',  # L227
  }  
  return parm

####################
## SKLearn machine learning models for DFG-motif conformation classification
def SKLearnDFGModelFiles():
  sk_dfg_model = {
    'rf': 'SK-0221_rf_model_dfg.pkl.bz2', 'svm': 'SK-0221_svm_rbf_model_dfg.pkl.bz2', 
    'nn': 'SK-0221_nn_model_dfg.pkl.bz2', 'kn':  'SK-0221_kn_model_dfg.pkl.bz2', 
    'gb': 'SK-0221_gb_model_dfg.pkl.bz2', 'gp':  'SK-0221_gp_model_dfg.pkl.bz2', 
    'dt': 'SK-0221_dt_model_dfg.pkl.bz2'  }
  return sk_dfg_model

## SKLearn machine learning models for C-helix/DFG conformation classification
def SKLearnKinfoModelFiles():
  sk_chx_model = {
    'rf': 'SK-0221_rf_model_full.pkl.bz2', 'svm': 'SK-0221_svm_lin_model_full.pkl.bz2', 
    'nn': 'SK-0221_nn_model_full.pkl.bz2', 'kn':  'SK-0221_kn_model_full.pkl.bz2', 
    'gb': 'SK-0221_gb_model_full.pkl.bz2', 'gp':  'SK-0221_gp_model_full.pkl.bz2', 
    'dt': 'SK-0221_dt_model_full.pkl.bz2'  }
  return sk_chx_model


##########################################################################
## Generate template setup file containing user-controlled variables
def GenerateTemplSetupScript( parameter_file, V=DefaultVariables() ):

  lst = [ 'SCRIPT', 'DATADIR','PDBDIR','HOMEDIR',None,
          'OUTPREF','PDBLIST','PDBALIGN',None,
          'USESKL',None,
          'CHECKPDB','CHKALIGN','SUPERPOSED',None,
          'MISSRES',None,'PYMOL','MPICPU' ]

  with open(parameter_file, 'w') as f:
    for key in lst:
      if key is None:
        f.write('\n')
      else:
        f.write('{0}\n{1}\t{2}\n\n'.format(V[key][1], key, V[key][0]))


##########################################################################
# Read in parameter file for residues and other params
def ParseParameterFile( param_list ):

  parm = DefaultVariables()

  for l in pd.read_csv(param_list, header=None, comment='#', sep='\s+').to_numpy():
    try:
      parm[l[0]][0] = l[1]
    except ValueError:
      sys.exit('\n \033[31m#2# Parameter Warning: Unknown parameter handle: \033[0m'+l[0])

  ## check for key items if they exist
  for flag in ['SCRIPT','DATADIR','PDBDIR','PDBLIST','PDBALIGN','PYMOL','REFPDB','MISSRES']:
    if flag == 'MISSRES':
      if re.search('none|false', parm['MISSRES'][0], re.IGNORECASE):
        parm['MISSRES'][0] = None
    elif not os.path.exists(parm[flag][0]):
      sys.exit('\n \033[31m#2# FATAL: Check \033[0m {0} - {1}'.format(flag,parm[flag][0]))

  ## adjust the on/off switch
  for flag in ['CHECKPDB','CHKALIGN','SUPERPOSED','USERRF']:
    if re.search('true', parm[flag][0], re.IGNORECASE):
      parm[flag][0] = True
    else:
      parm[flag][0] = False

  ## number of CPU for MPI
  try:
    parm['MPICPU'][0] = int(parm['MPICPU'][0])
  except TypeError:
    parm['MPICPU'][0] = 0
    print('\n \033[31mERROR: "MPICPU" must be integer - use all CPU instead: \033[35m{0}\033[0m'.format(parm['MPICPU'][0]))

  return parm


##########################################################################
