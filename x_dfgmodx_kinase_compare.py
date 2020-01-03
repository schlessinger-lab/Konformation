#!/usr/bin/env python3

import sys,re,os

import numpy as np
import pandas as pd

from Bio import SeqIO

from CommonUtility import file_handle

##########################################################################
# Calculate sequence identity and similarity of a query seq to a library of 
# sequence (or single seq) and output a list with the best one at the 1st row
### Original from --> 3_DFGmodx/x_homolog_templ_check.py
def BlastpPairwiseIdentity( result_directory, mdl_prot_fasta, kinase_profile ):

  # If input Fasta is a file, reconfigure to only the fasta name  
  if os.path.isfile(mdl_prot_fasta):
    fasta_name = mdl_prot_fasta.split('.fasta')[0]
  else:
    fasta_name = mdl_prot_fasta

#  print('\n  ** Calculate Sequence Identity between Query and Profile Sequences **')
#  print('  Query Fasta:    '+fasta_name+'.fasta')
#  print('  Kinase Profile: '+kinase_profile)
  # blastp to output: Name, AA_length, percent identity, percent positive
  # result in .csv format, omit other irrelevant data
#  print('blastp -query "{0}" -subject "{1}" -max_target_seqs 5000 -out {2}/{3}.idmat.txt -outfmt "6 sseqid length pident ppos"'.format(fasta_name+'.fasta', kinase_profile, result_directory, fasta_name.split('/')[-1]))
  os.system('blastp -query "{0}" -subject "{1}" -max_target_seqs 5000 -out {2}/{3}.idmat.txt -outfmt "6 sseqid length pident ppos"'.format(fasta_name+'.fasta', kinase_profile, result_directory, fasta_name.split('/')[-1]))

  # Parse percent identity result generated by BlastP. Did not use clustalo or
  # t_coffee because they do redundent pairwise identity calculation for other 
  # kinases to create a true matrix and that is not needed; only need 1 set of
  # pairwise identity between query sequence and the database sequences
  if not os.path.isfile('{0}/{1}.idmat.txt'.format(
                        result_directory, fasta_name.split('/')[-1])):
    print('\n  > \033[31m#2#\033[0m Alignment Warning: No Blastp output. Seq identity too low? '+fasta_name)
    return None
  elif os.stat('{0}/{1}.idmat.txt'.format(
                result_directory, fasta_name.split('/')[-1])).st_size == 0:
    print('\n  > \033[31m#2#\033[0m Alignment Warning: Blastp failed. Seq Identity too low? '+fasta_name)
    return None


  ## Extract the identity information from Blastp result; sometimes a chain is
  ## broken into fragments, need to combine them according to residue ratios
  Ident = {}
  with open('{0}/{1}.idmat.txt'.format(result_directory, fasta_name.split('/')[-1]), 'rU') as fi:
    for line in fi:
      Items = line.split('\t')
      name, aa, identity, positive = (Items[0].split('|')[0], int(Items[1]), 
                                      float(Items[2]), float(Items[3]) )
      if name in Ident:
        Ident[name].append([name, aa, identity, positive])
      else:
        Ident[name] = [ [name, aa, identity, positive] ]

  # Convert dictionary into Tulip data. If a Fasta name has multiple lines, 
  # the alignment/identity calculation is broken down into pieces for 1 seq.
  # Summarize the pieces into 1 by adding up the ratio
  Data = []
  for name in Ident:
    length = sum(list(zip(*Ident[name]))[1])  # rearrange tulip groups 
    x, y = 0.0, 0.0
    nm   = name.split('_')
    if enumerate(nm) != 2:
      nm.append('A')
    for row in Ident[name]:
      x += row[1] * row[2]
      y += row[1] * row[3]

    Data.append( [nm[0], nm[1], length, (x/length), (y/length)] )

  # sort the dataset by percent identity or positive, then by available length,
  # then by filename to prefer A or B, etc
  pdata = pd.DataFrame(Data)  
  pdata.columns = ['pdb_id', 'chain', 'length', 'identity', 'similarity']  
  pdata['pdb_full'] = pdata['pdb_id']+'_'+pdata['chain']

  pdata = pdata.sort_values( by=['identity', 'length', 'chain'], 
                              ascending=[False, False, True] )
  pdata_temp = pdata.drop('pdb_id',1).drop('chain',1)
  col = pdata_temp.columns.tolist()
  col = col[-1:] + col[:-1]
  pdata_temp = pdata_temp[col]
  pdata = pdata_temp

  pdata.to_csv('{0}/{1}.idmat.sort.txt'.format(result_directory, fasta_name.split('/')[-1]), sep='\t', encoding='utf-8', float_format='%4.2f', index=False)

  Data = [[r.pdb_full, r.length, r.identity, r.similarity ] for idx, r in pdata.iterrows()]

  return Data


##########################################################################
