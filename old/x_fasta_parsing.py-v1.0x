#!/usr/bin/python

import re,os,glob,sys
import numpy as np
from CommonUtility import *
from aa_residue import AA
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
p = PDBParser(PERMISSIVE=1)

##########################################################################
## Check sequence if there is empty space ' ' or '-'
def CheckSequence(Seq):
  for res in Seq:
    if res == '-' or res == ' ':
      return False
  return True

##########################################################################
## Get the sequence from the reference PDB file
def RefFastaColumn(Fasta_Lib, residues, name):

  seq, Set = '', {}
  for seq_record in Fasta_Lib:
    if re.search(r'{0}'.format(name), seq_record.id.replace(':','_')):
      seq = str(seq_record.seq)

  # Find the sequence by the 5/7-resid code in (D-2,D-1,D,D+1,D+2) format
  # but refer the column by the center residue 'D'
  Match  = re.finditer(r'{0}'.format(residues), seq)
  for m in Match:
    column = m.start(0)
  Set[name.split('.')[0]] = list(residues)
  print(' reference sequence:\t'+residues)
  print(' Column in Fasta:\t'+str(column+1))  # convert from 'start 0' to 1
  return Set, column


##########################################################################
## Read from alignment FASTA. Find the position of the search pattern in the 
## sequence and make a new sequence, with the residue before the input sequence
## plus the first two residues in the input sequence.
## ** Make sure no 'blank' residue in fasta library in that region
def FastaFromColumn(Fasta_Lib, column, res_num):
  Set    = {}
  for seq_record in Fasta_Lib:
    # Since all query sequences are aligned, the column numbers should 
    # correspond to the correct helix in query sequences
    seq    = str(seq_record.seq)
    Qu_Seq = [seq[column+i] for i in range(0, res_num)]
    pdb_id = seq_record.id.split('|')[0]
    Set[pdb_id.replace(':','_')] = Qu_Seq
  return Set


#########################################################################
#
#   Peter M.U. Ung  @ MSSM
#   
#   v0.1    17.01.09
#   v0.2    17.01.28
#   v0.3    17.01.30
#
#   v1.0    17.02.02
#
