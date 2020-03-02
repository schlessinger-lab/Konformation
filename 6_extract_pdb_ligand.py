#!/usr/bin/env python3

import sys,os

##########################################################################
#
#  Peter M.U. Ung @ MSSM/Yale
#
#  v1.0  20.03.01
#
#  From an input file of Parsed PDB Header file with information on 
#  "pdb_id" and "ligand", download the SMILES string of bound ligands
#  and associate these ligands to the "pdb_id" they came from.
#
##########################################################################
msg = '''
  > {0}
\t[ Input file of Parsed PBB Header file with "pdb_id" and "ligand" ]
\t[ Output prefix ]
'''.format(sys.argv[0])
if len(sys.argv) != 3: sys.exit(msg)

import re
import time
import requests
import xmltodict

import numpy as np
import pandas as pd

from rdkit import Chem

##########################################################################
def main( infile, outpref ):

  ## Read in kinase data
  pdb_df = pd.read_csv(infile, comment='#', sep=',')
  df = pdb_df[['pdb_id','chain_id','ligand']].dropna().to_numpy()

  lig_dict = [ DownloadLigand(lig) for lig in RearrangeKinaseData(df).items() ]
  columns  = ['lig_id','pdb_id','formula','lig_mw','lig_name','smiles']
  lig_df   = pd.DataFrame(lig_dict, columns=columns)

  lig_df.to_csv(outpref+'.csv.gz', sep=',',index=None)
  lig_df.to_excel(outpref+'.xlsx', index=None)


##########################################################################
## Parse kinase data to get unique ligand PDB_ID and associated PDB(s)
def RearrangeKinaseData( df ):
  lst = []

  for row in df:
    if row[2] is np.nan:
      continue
    if re.search('|', row[2]):
      for s in row[2].split('|'):
        tmp = {}
        tmp['pdb_id'] = row[0]+'_'+row[1]
        tmp['lig_id'] = s.split(':')[0]
        lst.append(tmp)
    else:
      tmp = {}
      tmp['pdb_id'] = row[0]+'_'+row[1]
      tmp['lig_id'] = row[2].split(':')[0]
      lst.append(tmp)

  new = {}
  for row in lst:
    if row['lig_id'] not in new:
      new[row['lig_id']] = ''.join(row['pdb_id'])
    else:
      if not re.search(row['pdb_id'], new[row['lig_id']], re.IGNORECASE):
        new[row['lig_id']] += '|'+row['pdb_id']

  return new


##########################################################################
## Parse Ligand data retreived from PDB_ID
def DownloadLigand( inp ):
  Rst = {'lig_id': inp[0], 'pdb_id': inp[1]}

  ## use 'requests' to download ligand
  lig_html = 'https://www.rcsb.org/pdb/rest/describeHet?chemicalID={0}'.format(Rst['lig_id'])
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

  return Rst

##########################################################################
if __name__ == '__main__':
  main( sys.argv[1], sys.argv[2] )
