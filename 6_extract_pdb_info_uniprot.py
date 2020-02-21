#!/usr/bin/python

##########################################################################
#
#   Peter M.U. Ung @ MSSM
#
#   v1.0 -  10.23.2014
#   v2.0    17.05.17    use MPI to run web search
#   v2.1    17.05.31    look up new info
#   v3.0    20.01.09    update data retrieve method for RCSB and Uniprot
#
#   Previous:   0_pdb_uniprot_length.py
#
#   Purpose:    From a list of PDB file (with chain ID associated in the name)
#               search the PDB chain length and the corresponding full-length
#               from RCSB PDB and UniProt databases via internet. No need to
#               read from hard-copy of PDB files.
#               e.g. PDB file name: 1ATP_A.xxx.pdb
#                   (PDB ID and chain ID separated by '_')
#
##########################################################################

import sys,os
msg = '''\n    {0}\n\t\t[list of PDB]\n\t\t[output prefix]\n'''.format(sys.argv[0])
if len(sys.argv) != 3: sys.exit(msg)

import re,glob
import pandas as pd

from x_extract_pdb_info import ProcessPDBList

##########################################################################
def main( pdb_list, output_pref ):

  print('  ## Search RCSB PDB and Uniprot webpages ##')
  Entity = ProcessPDBList(pdb_list, output_pref)

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
#        print(dict(Prot[key]))
        if re.search(Prot[key].chain_id, key):
          Rst.append(dict(Prot[key]))
        continue

  print('  ## Print out data: {0} ##'.format(len(Rst)))  
  df = pd.DataFrame(Rst)
  df.index.name = 'index'

  columns = [ 'pdb_id','chain_id','chain','pdb_length','uni_length',
              'uni_id','gene','p_name','mutate','mutation','ec',
              'species','common','taxid','deposit','release','latest',
              'ligand','aa_modif','resolu','space','pmid' ]

  reorder = df[columns]

  # Output to Excel
  reorder.to_excel(output_pref+'.xlsx', header=True , na_rep='NaN')

  # Output to CSV
  reorder.to_csv( output_pref+'.csv.gz', sep=',', na_rep='NaN', header=True,
                  encoding='utf-8', ) #cols=columns)

  # Output to Pickle
  df.to_pickle(output_pref+'.pkl', compression='bz2')

  os.system("rm _TEMP*.HTML")


##########################################################################
##########################################################################
if __name__ == "__main__":
  main(sys.argv[1], sys.argv[2])
