#!/usr/bin/env python3

import sys,os,re
import pandas as pd

def DefaultVariables():
  script   = '/home/pmung/Dropbox/9_scripts/3_program/structures/4_Konformation/'
  data_dir = script+'z_database/'
  outpref  = 'check_kinase.NEWEST'

  Settings = {
    'WORKDIR':   [os.getcwd(), '# Current working directory'],
    'OUTPREF':   [outpref, '# Output prefix for result and intermediate info files'],

    'RDKNOWNKIN':[data_dir+'all_downloaded_kinase_pdb.200220.list','# Read from list of downloaded kinase structures (with extension .1atp.pdb)'],
    'RDNONKIN':  ['check_kinase.200220.non_kinase.list','# Read from a List of PDBs that has no kinase catalytic domain'],
    'RDCHECKED': ['check_kinase.200220.checked_pdb.list', '# Read from a List of any PDB that have been checked previously to skip redundant checks'],

    'WTNEWKNOWN':['all_downloaded_kinase_pdb.NEWEST.list','# Write to List of newly found kinase structures, including KNOWNLIST'],
    'WTNONKIN':  [outpref+'.non_kinase.list','# Write to a List of PDBs that has no kinase catalytic domain'],
    'WTCHECKED': [outpref+'.checked_pdb.list', '# Write to a List of PDB that have been checked for kinase domain'],

    'KINOMEGPDB':[data_dir+'MD_human_kinome_alignment.2019-2.fasta','# FASTA file of aligned canonical human kinome, with gap'],
    'KINOMEDB':  [data_dir+'MD_human_kinome_alignment.2019.nogap.fasta','# FASTA file of unaligned canonical human kinome, NO gap'],
    
    'LENCUTOFF': [225,'# Sequence length Cutoff for kinase catalytic domain (def: 220)'],
    'IDTCUTOFF': [40.,'# Sequence identity cutoff for kinase catalytic domain recognition (def: 40.0%)'],

    'REFPDB':    [data_dir+'1ATP_E.pdb','# Reference PDB, bovine PKA (1ATP_E)'],
    'REFRES':    ['"resi 121-139+162-183"','# Reference 1ATP_E residues for PyMOL superposition'],
    'OUTEXT':    ['1atp.pdb','# Output extension for 1atp-superposed structure'],
    'PYMOL':     ['/usr/bin/pymol','# Full path to PyMOL Executable'],
  }
  return Settings


##########################################################################
## Write out a template setting file for read-in
def GenerateTemplSetupScript( parameter_file, Vars=DefaultVariables() ):
  with open(parameter_file, 'w') as fo:
    for key in list(Vars.keys()):
      fo.write('{0}\n{1}\t{2}\n\n'.format(Vars[key][1], key, Vars[key][0]))


##########################################################################
# Read in parameter file for residues and other params
def ParseParameterFile( param_list ):

  Settings = DefaultVariables()
  for key in Settings:
    Settings[key] = Settings[key].pop(0)

  items = pd.read_csv(param_list, comment='#', header=None, sep='\s+')

  for idx, row in items.iterrows():
    if row[0] in Settings:
      Settings[row[0]] = row[1]
    else:
      print('\n  \033[31m#2#\033[0m Parameter Warning: Unknown parameter handle: \033[31m{}\033[0m'.format(row[0]))
  return Settings


##########################################################################
#
#  This script provides the Default Setting for 5_update_kinase_db.py
#
##########################################################################