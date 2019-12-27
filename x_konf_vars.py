#!/usr/bin/python

import sys,os,re

def DefaultVariables():
  pdb_dir = '/home/pmung/Dropbox/1_kinase/1_family/1_stdkinases/170109/'
  script  = '/home/pmung/Dropbox/9_scripts/3_program/structures/4_Konformation/'
  data_dir= script+'z_database/'

  RefRes = {
    'HOMEDIR':  os.getcwd(),
    'PDBLIST':  None,   # List of input PDB structure to examine
    'SUPERPO':  'False',  # Is input PDB pre-superposed to 1ATP ref residue
    'MISSRES':  None,   # User-input list of missing residue PDB

    'PDBDIR':   pdb_dir+'2_align', # Directory to internal PDB for referencing
    'MISSREF':  pdb_dir+'3_kinase_param/stdy_kinase.param.170825.missing.root.txt',  # List of residue PDB known to missed out in the 5-resid-based search algorithm
    'FASTA':    data_dir+'stdy_kinase.clean_3.180709.tcoffee_1d.fasta', # Manually corrected FASTA MSA of all internal PDB for referencing
    'BLASTDB':  data_dir+'stdy_kinase.clean_3.180709.nogap.tcoffee_1d.blastdb', # BlastDB version of FASTA of all internal PDB, with no alignment/gaps
    'SCRIPT':   script, # Directory of all relevant scripts

    'REFPDB':   data_dir+'1ATP_E.pdb',   # Reference PDB, bovine PKA (1ATP_E)
    'REFRES':   'resi 122-138+162-183', # reference 1ATP_E PDB residues for PyMOL superposition
    'REMOVE':   None,
    'OUTEXT':   '.1atp.pdb',   # Output extension for all 1atp-superposed PDB
  
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
    f.write('{0}\t\t{1}\n'.format('PDBLIST', '*.list'))
    f.write('{0}\t\t{1}\n'.format('SUPERPO', V['SUPERPO']))
    f.write('{0}\t\t{1}\n'.format('MISSRES', V['MISSRES']))


##########################################################################
# Read in parameter file for residues and other params
def ParseParameterFile( param_list ):

  RefRes = DefaultVariables()

  with open(param_list, 'rh') as fi:
    lines = filter(None, (l.rstrip() for l in fi))
  for l in lines:
    if re.search(r'^#', l): continue
    i = l.split()
    try:
      RefRes[i[0]] = i[1]
    except ValueError:
      print('\n  #2# Parameter Warning: Unknown parameter handle: '+i[0])
      continue
  
  if RefRes['PDBLIST'] is None:
    sys.exit('\n  #2# Parameter FATAL: Input PDB list is not available ##')
  if type(RefRes['MISSRES']) is str: 
    if re.search(r'None', RefRes['MISSRES'], re.IGNORECASE):
      RefRes['MISSRES'] = None
  if type(RefRes['SUPERPO']) is str:
    if re.search(r'true', RefRes['SUPERPO'], re.IGNORECASE):
      RefRes['SUPERPO'] = True
    else:
      RefRes['SUPERPO'] = False

  return RefRes


##########################################################################
