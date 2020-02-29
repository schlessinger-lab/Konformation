#!/usr/bin/env python3

##########################################################################
#
#   Peter Man-Un Ung @ MSSM
#
#   v1.0    17.10.07	count how many kinase in each conformational state
#
#   Read in the classification results after R machine learning and parse
#   to count how many unique kinase (uniprot_id) in each conformation state
#
#   format: gene_name, species, unip_id, conf_num, conf, pdb_num, pdb_list
#
##########################################################################

import pandas
import sys, os, re

msg = ''' '''
if len(sys.argv) != 2: sys.exit(msg)

##########################################################################

## read database(.csv) into dataframe
Data = pandas.read_csv(sys.argv[1], sep=',', header='infer', )

Result = {}
group  = {}

for idx, row in Data.iterrows():
  group = {}

  if len(str(row['gene']).split()) > 1:
    line = re.sub(r',', '', str(row['gene']))
    gene = line.split()[0]
  else:
    gene = str(row['gene'])

  species = re.sub(r' ', '_', str(row['species']))
  
  if row['uni_id'] in Result:
    if row['group'] in Result[row['uni_id']]:
      Result[row['uni_id']][row['group']][0].append(row['pdb_id'])
    else:
      Result[row['uni_id']][row['group']] = [[row['pdb_id']], gene, species]

    if not re.search(r'nan|NaN|NAN', gene):
      for idx, group in enumerate(Result[row['uni_id']]):
        Result[row['uni_id']][group][1] = gene

  else:
    group[row['group']] = [[row['pdb_id']], gene, species]
    Result[row['uni_id']] = group


for unip in Result:
  grp_num = len(Result[unip])
  for group in Result[unip]:
    pdb_num = len(Result[unip][group][0])
    pdb_lst = ' '.join( Result[unip][group][0] )
    gene = Result[unip][group][1]
    spec = Result[unip][group][2]
    if grp_num > 0:
      print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format(
              gene, spec, unip, grp_num, group, pdb_num, pdb_lst ))
