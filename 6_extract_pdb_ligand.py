#!/usr/bin/env python3

import sys,os

##########################################################################
#
#  Peter M.U. Ung @ MSSM/Yale
#
#  v1.0  20.03.01
#
#  From input file(s) of Parsed PDB Header file with information on 
#  "pdb_id" and "ligand", and/pr "conf" information (e.g. result from
#  ) download the SMILES string of bound ligands
#  and associate these ligands to the "pdb_id" they came from.
#
##########################################################################
msg = '''
  > {0}
\t-in <file1> <file2>  [ File(s) of Parsed PBB Header file with "pdb_id" and "ligand", and/or "conf" ]
\t-op <str>            [ Output prefix ]
'''.format(sys.argv[0])
if len(sys.argv) == 1: sys.exit(msg)

import re
import time
import requests
import xmltodict

import numpy as np
import pandas as pd

from rdkit import Chem
from collections import Counter
from argparse import ArgumentParser

##########################################################################
def main( ):
  args  = cmdlineparse()

  ## Read in multiple kinase data
  files  = [ pd.read_csv(inf, comment='#', sep=',') for inf in args.infiles ]
  pdb_df = pd.concat(files, ignore_index=True)
  if 'pdb' not in pdb_df.columns:
    pdb_df['pdb'] = np.nan
  if 'chain_id' not in pdb_df.columns:
    pdb_df['chain_id'] = np.nan
  if 'conf' not in pdb_df.columns:
    pdb_df['conf'] = np.nan
  df = pdb_df[['pdb_id','pdb','chain_id','conf','ligand']].dropna().to_numpy()

  Ligs = RearrangeKinaseData(df)
  print('  # Number of Ligand found: \033[31m{0}\033[0m'.format(len(Ligs)))

  columns  = ['lig_id','pdb_id','cidi','cido','codi','codo','wcd','formula','lig_mw','lig_name','smiles']
  lig_dict = [ DownloadLigand(lig, lst, columns) for lig, lst in Ligs.items() ]
  lig_df   = pd.DataFrame(lig_dict, columns=columns)

  lig_df.to_csv(args.outpref+'.csv.gz', sep=',',index=None)
  lig_df.to_excel(args.outpref+'.xlsx', index=None)


##########################################################################
## Parse Ligand data retreived from PDB_ID
def DownloadLigand( lig, lst, columns ):

  Rst = { 'lig_id': lig, 'pdb_id': lst[0], 
          'cidi':0, 'cido':0, 'codi':0, 'codo':0, 'wcd':0 }

  y = Counter(lst[1].split('|'))
  for conf in y.keys():  Rst[conf] = y[conf]

  ## use 'requests' to download ligand
  lig_html = 'https://www.rcsb.org/pdb/rest/describeHet?chemicalID={0}'.format(lig)
  lig_data = requests.get(lig_html)

#  time.sleep(0.1)

  url_dict = xmltodict.parse(lig_data.content.decode(), process_namespaces=True)
  ligx = url_dict['describeHet']['ligandInfo']['ligand']

  Rst['formula']  = ligx['formula']
  Rst['lig_mw']   = ligx['@molecularWeight']
  Rst['lig_name'] = ligx['chemicalName']
  try:
    Rst['smiles'] = Chem.MolToSmiles(Chem.MolFromSmiles(ligx['smiles']))
  except TypeError:
    Rst['smiles'] = ligx['smiles']

  Series = [ Rst[name] for name in columns ]
  return Series


##########################################################################
## Parse kinase data to get unique ligand PDB_ID and associated PDB(s)
def RearrangeKinaseData( df ):
  lst = []

  # convert to list of dict
  # ['pdb_id','pdb','chain_id','conf','ligand']
  for inf in df:
    if inf[4] is np.nan:
      continue
    if re.search('|', inf[4]):
      for s in inf[4].split('|'):
        row = {}
        row['pdb_id'] = inf[0]
        row['conf']   = inf[3]
        row['lig_id'] = s.split(':')[0]
        lst.append(row)
    else:
      row = {}
      row['pdb_id'] = inf[0]
      row['conf']   = inf[3]
      row['lig_id'] = inf[4].split(':')[0]
      lst.append(row)

  # convert list of data to dict of list
  Ligs = {}
  for row in lst:
    if row['lig_id'] not in Ligs:
      pdb_id = ''.join(row['pdb_id'])
      conf   = ''.join(row['conf'])
      Ligs[row['lig_id']] = [ pdb_id, conf ]
    else:
      if not re.search(row['pdb_id'], Ligs[row['lig_id']][0], re.IGNORECASE):
        pdb_id = Ligs[row['lig_id']][0]+'|'+row['pdb_id']
        conf   = Ligs[row['lig_id']][1]+'|'+row['conf']
        Ligs[row['lig_id']] = [ pdb_id, conf ]
      else:
        pdb_id = Ligs[row['lig_id']][0]
        conf   = Ligs[row['lig_id']][1]+'|'+row['conf']
        Ligs[row['lig_id']] = [ pdb_id, conf ]

  return Ligs


##########################################################################
def cmdlineparse():
  p = ArgumentParser(description="command line arguments")

  p.add_argument('-in', dest='infiles', required=True, nargs='+',
                      help='One or more file of Parsed PBB Header file with "pdb_id" and "ligand"')
  p.add_argument('-op', dest='outpref', required=True,
                  help='Output prefix')

  return p.parse_args()

##########################################################################
if __name__ == '__main__':
  main(  )
