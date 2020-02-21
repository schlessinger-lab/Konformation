#!/usr/bin/env python3

import re
import os

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.PDB.PDBParser   import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

from aa_residue import AA
from x_variables import per_line
from x_pir_edit  import CheckPIR
from x_pir_edit  import ModifyPIR

##############################################################################

##########################################################################
## Alignment using a pre-existing MSA profile. T-Coffee has issue with this. 
## MUSCLE and ClustalO came out about same time, 2004 and 2003, MUSCLE performs
## the best especially when doing single-seq profile alignment. Output the 
## profile-aligned fasta object
def MuscleProfileAlign( fasta_database, fasta_file, temp_file ):

  os.system('muscle -profile -in1 "{0}" -in2 "{1}" -out "{2}" -maxiters 64 -quiet'.format(
              fasta_database, fasta_file, temp_file ) )

  Tget_List = list(SeqIO.parse(temp_file, 'fasta'))

  return Tget_List[1]


##############################################################################
# Reformat and no alignment: Remove empty gap column from aligned FASTA file  
def RemoveFastaGapColumn( fasta_input, fasta_output ):

  os.system('t_coffee -other_pg seq_reformat -in "{0}" -action +rm_gap -output=fasta > {1}'.format(fasta_input, fasta_output))


##########################################################################
