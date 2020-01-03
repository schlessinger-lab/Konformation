#!/usr/bin/env python3

import sys,os,re
import pandas as pd

def DefaultVariables():
  script   = '/home/pmung/Dropbox/9_scripts/3_program/structures/4_Konformation/'
  data_dir = script+'z_database/'

  Settings = {
    'WORKDIR':   [os.getcwd(), '# Current working directory'],
    'OUTPREF':   ['test', '# Output prefix for intermediate info files'],

    'KNOWNLIST': ['downloaded_stdy_kinase.170801.list','# Read from list of downloaded kinase structures (*.1atp.pdb)'],
    'NEWKNOWN':  ['downloaded_stdy_kinase.191230.list','# Save to updated kinase structure list'],

    'KINASEDB':  [data_dir+'human_kinome_alignment.nogap.pmung.180614.fasta','# FASTA file of human kinome'],
#    'NONKINASE': [data_dir+'non_kinase_pdb.191231.list','# List of PDBs that has no kinase catalytic domain'],
    'NONKINASE': ['None', '###'],
    'CHECKEDPDB':['None', '###'],

    'LENCUTOFF': [220,'# Sequence length Cutoff for kinase catalytic domain (def: 220)'],
    'IDTCUTOFF': [40.,'# Sequence identity cutoff for kinase catalytic domain recognition (def: 40.0%)'],

    'REFPDB':    [data_dir+'1ATP_E.pdb','# Reference PDB, bovine PKA (1ATP_E)'],
    'REFRES':    ['"resi 122-138+162-183"','# Reference 1ATP_E residues for PyMOL superposition'],
    'OUTEXT':    ['.1atp.pdb','# Output extension for 1atp-superposed structure']
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