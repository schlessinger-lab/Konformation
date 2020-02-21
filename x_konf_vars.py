#!/usr/bin/env python3

import sys,os,re
import pandas as pd

def DefaultVariables():
  pdb_dir = '/home/pmung/Dropbox/1_kinase/1_family/1_stdkinases/170109/'
  script  = '/home/pmung/Dropbox/9_scripts/3_program/structures/4_Konformation/'
  data_dir= script+'z_database/'

  RefRes = {
    'PDBDIR':   pdb_dir,     # Path to directory with stored PDB structures
    'SCRIPT':   script,      # Path to main scripts
    'DATADIR':  data_dir,    # Path to main scripts' database
    'HOMEDIR':  os.getcwd(),    # Path to home directory for working

    'PDBLIST':  'pdb.list',     # List of input PDB structure to examine
    'CHKPDB':   'True',   # Check input PDB provided in PDBLIST
    'CHKALIGN': 'True',   # Check alignment of input PDB seq to remove non-kinases
    'SUPERPOSED': 'False',  # Has input PDB pre-superposed to 1ATP ref residue
    'MISSRES':  None,     # User-input list of missing residue PDB

    'PDBDIR':   pdb_dir+'2_align', # Directory to internal PDB for referencing
    'MISSREF':  pdb_dir+'3_kinase_param/stdy_kinase.param.170825.missing.root.txt',  # List of residue PDB known to missed out in the 5-resid-based search algorithm
    'FASTA':    data_dir+'stdy_kinase.clean_3.180709.tcoffee_1d.fasta', # Manually corrected FASTA MSA of all internal PDB for referencing
    'BLASTDB':  data_dir+'stdy_kinase.clean_3.180709.nogap.tcoffee_1d.blastdb', # BlastDB version of FASTA of all internal PDB, with no alignment/gaps
    'SCRIPT':   script, # Directory of all relevant scripts

    'REFPDB':   data_dir+'1ATP_E.pdb',   # Reference PDB, bovine PKA (1ATP_E)
    'REFRES':   'resi 121-139+162-183',  # Reference 1ATP_E PDB residues for PyMOL superposition
    'OUTEXT':   '.1atp.pdb',   # Output extension for all 1atp-superposed PDB
    'OUTPREF':  'output_prefix',    # Output prefix for data generated

## 5-resid based reference sequences for searching other kinase PDB's relevant
## residue positions. Intended residue is always at the center of the sequence
## flanked by equal no. of residue on both sides. Need at least 5-resid to have
## good confidence in correct identification. Search algorithm will fail if PDB
## structure has missing or unnatural residue in corresponding sequence. Will
## need to manually find and save the missing residue into a MISSING PDB, and 
## append the info into the MISSING list. 
## For failed cases, try alternative 5-resid sequence that takes last residue

    'HELIX':    'HTLNEKRIL',  # [Ref C-helix (Glu) | 5,7,9-resid center on E91]
    'NDOM':     'AMKIL',  # [Ref N-lobe beta3(Lys)] K72
    'CDOM':     'VTDFG',  # [Ref C-lobe DFG(Asp)] D184
    'DFGF':     'TDFGF',  # [Ref C-lobe DFG(Phe)] F185

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
  return RefRes

##########################################################################
def GenerateTemplSetupScript( parameter_file, V=DefaultVariables() ):

  with open(parameter_file, 'w') as f:
    f.write('{0}\t\t{1}\n'.format('SCRIPT', V['SCRIPT']))
    f.write('{0}\t\t{1}\n'.format('DATADIR', V['DATADIR']))
    f.write('{0}\t\t{1}\n'.format('PDBDIR', V['PDBDIR']))
    f.write('\n')
    f.write('{0}\t\t{1}\n'.format('PDBLIST', '*.list'))
    f.write('\n')
    f.write('{0}\t\t{1}\n'.format('CHKPDB', V['CHKPDB']))
    f.write('{0}\t\t{1}\n'.format('CHKALIGN', V['CHKALIGN']))
    f.write('{0}\t\t{1}\n'.format('SUPERPOSED', V['SUPERPOSED']))
    f.write('{0}\t\t{1}\n'.format('MISSRES', V['MISSRES']))
    f.write('\n')
    f.write('{0}\t\t{1}\n'.format('OUTPREF', V['OUTPREF']))


##########################################################################
# Read in parameter file for residues and other params
def ParseParameterFile( param_list ):

  RefRes = DefaultVariables()

  lines = pd.read_csv(param_list, header=None, comment='#', sep='\s+').values.tolist()
  for l in lines:
    try:
      RefRes[l[0]] = l[1]
    except ValueError:
      sys.exit('\n  #2# Parameter Warning: Unknown parameter handle: \033[31m{}\033[0m'.format(l[0]))
      continue

  if not os.path.exists(RefRes['SCRIPT']):
    sys.exit('\n  #2# FATAL: Directory is not available: \033[31mSCRIPT\33[0m')
  if not os.path.exists(RefRes['DATADIR']):
    sys.exit('\n  #2# FATAL: Directory is not available: \033[31mDATADIR\033[0m')
  if not os.path.exists(RefRes['PDBDIR']):
    sys.exit('\n  #2# FATAL: Directory is not available: \033[31mPDBDIR\033[0m')

  if RefRes['PDBLIST'] is None:
    sys.exit('\n  #2# Parameter FATAL: Input is not available: \033[31mPDBLIST\033[0m')
  if RefRes['OUTPREF'] is None:
    sys.exit('\n  #2# Parameter FATAL: Input is not available: \033[31mOUTPREF\033[0m')
  
  if type(RefRes['MISSRES']) is str: 
    if re.search('none', RefRes['MISSRES'], re.IGNORECASE):
      RefRes['MISSRES'] = None
  if type(RefRes['CHKALIGN']) is str:
    if re.search('none', RefRes['CHKALIGN'], re.IGNORECASE):
      RefRes['CHKALIGN'] = None

  if type(RefRes['SUPERPOSED']) is str:
    if re.search('true', RefRes['SUPERPOSED'], re.IGNORECASE):
      RefRes['SUPERPOSED'] = True
    else:
      RefRes['SUPERPOSED'] = False
  if type(RefRes['CHKPDB']) is str:
    if re.search('true', RefRes['CHKPDB'], re.IGNORECASE):
      RefRes['CHKPDB'] = True
    else:
      RefRes['CHKPDB'] = False


  return RefRes


##########################################################################
