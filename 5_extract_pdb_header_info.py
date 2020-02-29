#!/usr/bin/env python3

##########################################################################
#
#   Peter M.U. Ung @ MSSM
#
#   v1.0 -  10.23.2014
#   v2.0    17.05.17    use MPI to run web search
#   v2.1    17.05.31    look up new info
#   v3.0    20.01.09    update data retrieve method for RCSB and Uniprot
#   v3.1    20.02.27    save salt ions and additives
#
#   Previous:   0_pdb_uniprot_length.py
#
#   From a list of PDB files (with chain ID associated in the name)
#   to extract structure information, ligand name, etc.
#   Search the PDB chain length and the corresponding full-length
#   from RCSB PDB and UniProt databases via internet. No need to
#   read from hard-copy of PDB files.
#
#   The search process cannot use MPI, RCSB blocks multiple requests
#   from same IP address simultaneous? (prevent DoS attack?)
#
#   Format of the list:
#     1) <pdb_id>_<chain_id>.xxx.pdb
#     2) <pdb_id>_<chain_id>.xxx.pdb <chain_id>
#     3) <pdb_id>.xxx.pdb            <chain_id>
#     4) <pdb_id>                    <chain_id>
#   (PDB ID and chain ID separated by '_')
#
#   e.g.:   1ATP_E.xxx.pdb
#           3HHP.xxx.pdb    C
#           6GTT            A
#
#                   
##########################################################################

import sys,os
msg = '''\n    {0}
\t\t[list of PDB]
\t\t[output prefix]\n
Format of the list:
  1) <pdb_id>_<chain_id>.xxx.pdb
  2) <pdb_id>_<chain_id>.xxx.pdb <chain_id>
  3) <pdb_id>.xxx.pdb            <chain_id>
  4) <pdb_id>                    <chain_id>
  (PDB ID and chain ID separated by '_')\n
  e.g.:   1ATP_E.xxx.pdb
          3HHP.xxx.pdb    C
          6GTT            A\n'''.format(sys.argv[0])
if len(sys.argv) != 3: sys.exit(msg)

import re,glob
import pandas as pd

from tqdm import tqdm

from x_extract_pdb_info import SearchPDB
from x_extract_pdb_info import ProcessPDBList

##########################################################################
def main( pdb_list, output_pref ):

  ## Read in the pdb list
  PDB = ProcessPDBList(pdb_list)
  print('Found entries: \033[31m{0}\033[0m'.format(len(PDB)))

  print('  ## Search RCSB PDB and Uniprot webpages ##')  
  ## The search process cannot use MPI, RCSB blocks multiple requests
  ## from same IP address simultaneous? (prevent DoS attack?)
  ## Followed by Extracting the UniProt ID of the PDB via internet
  print(PDB)
  Entity = [SearchPDB(pdb) for pdb in tqdm(PDB)]
  print('\n  -- Search Completed for \033[31m{0}\033[0m / \033[31m{1}\033[0m entries\n'.format(len(PDB), len(Entity)))


  ## Reformat the Entity data object to dictionary for dataframe construction
  Rst = []
  for Prot in Entity:
    for key in Prot.keys():
      if Prot[key].chain_id is None or Prot[key].chain is None:
        print('  \033[31mWarning:\033[0m Found PDB with missing chain: \033[31m{}\033[0m'.format(key))
        print(dict(Prot[key]))
        print(Prot[key].chain_id)
        print(Prot[key].chain)
        continue
      else:
        if re.search(Prot[key].chain_id, key):
          Rst.append(dict(Prot[key]))
        continue


  print('  ## Print out data: \033[31m{0}\033[0m ##'.format(len(Rst)))  
  df = pd.DataFrame(Rst)
  df.index.name = 'index'

  columns = [ 'pdb_id','chain_id','pdb_length','uni_length',
              'uni_id','gene','p_name','mutate','mutation','ec',
              'species','common','taxid','deposit','release','latest',
              'ligand','salt','aa_modif','resolu','space','pmid' ]

  reorder = df[columns]

  ## Output to Excel
  reorder.to_excel(output_pref+'.xlsx', header=True , na_rep='NaN', index=None)

  ## Output to CSV
  reorder.to_csv( output_pref+'.csv.gz', sep=',', na_rep='NaN', header=True,
                  index=None, encoding='utf-8' )

  ## Output to Pickle
#  reorder.to_pickle(output_pref+'.pkl.bz2', compression='bz2')


##########################################################################
##########################################################################
if __name__ == "__main__":
  main(sys.argv[1], sys.argv[2])
