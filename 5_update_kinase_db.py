#!/usr/bin/env python3

import sys,os

##########################################################################
#
#  Peter MU Ung @ MSSM/Yale
#
#  v1  19.12.28
#  v2  20.01.12  generate MSA alignment for all new PDB for Konformation
#  v3  20.02.20  clean up read/write of result files
#
#  Search RCSB Protein Data Bank for all kinase-associated structures with
#  Queries, then compare the new query results to existing internal kinase
#  1ATP-superposed PDB library (known_list) to find those that aren't in the
#  library yet. These new entries are checked (seq length, seq identity to
#  human kinome) before downloading the kinase structures, which are then
#  superposed onto 1ATP using a standard reference residue set. The new 
#  results are added to the existing internal kinase library list (known_list)
#
#  * as long as the (known_list) is kept current and the downloaded PDBs
#    are placed into the correct folder, this script can be made into a
#    repeating "cron" job to update the kinase structure library.
#
#  * tried using MPI to speed up the download of FASTA and PDB but seems to
#    run into server dockdown if multiple jobs are retrieving from the same
#    IP. Could be a defense mechanism against DOS attack?? Using serial
#    retrieve with a delay of 0.1s in between seem to work fine.
#
##########################################################################
import re
import io
import time
import requests
import contextlib

import pandas as pd

from Bio import PDB
from Bio import SeqIO
from tqdm import tqdm
from pathos import multiprocessing
from argparse import ArgumentParser

from x_pdb_superpose import SuperposePDB
from x_query_var import ParseParameterFile
from x_search_align import CacheSeqDatabase
from x_search_align import RemoveFastaGapColumn
from x_query_var import GenerateTemplSetupScript
from x_search_align import GenerateProfileAlignment
from x_kinase_pdb_search_query import DownloadFASTA
from x_kinase_pdb_search_query import DownloadNewPDB
from x_kinase_pdb_search_query import CheckExistingPDBs
from x_kinase_pdb_search_query import ReadKnownKinaseList
from x_kinase_pdb_search_query import SearchRCSBWithQuery
from x_kinase_pdb_search_query import KinaseDownloadQueries
from x_kinase_pdb_search_query import CheckKinaseSeqIdentity


#######################################################################
msg = '''
  > {0}
\t-r <file>  [ Read in Parameter file ]\n
\t-set       [ Optional: Generate Parameter file used by "-r" ]
\t-p <file>  [ Optional: List of PDB (and Chain ID) to directly download (format: pdb_id\tchain_id) ]
'''.format(sys.argv[0])
if len(sys.argv) == 1: sys.exit(msg)

def main():
  args  = cmdlineparse()

  if args.set:
    if args.read_param:
      GenerateTemplSetupScript(args.read_param)
      sys.exit('\n  \033[34m## Printing out Setting File. Rerun with -r <setting file> to run ##\033[0m')
    else:
      GenerateTemplSetupScript('kinase_query.template.setup')
      sys.exit('\n  \033[34m## Printing out Setting File. Rerun with -r <setting file> to run ##\033[0m')
  if args.read_param:
    parms = ParseParameterFile(args.read_param)
    print('\n  \033[34m## Reading in Parameter file ##\033[0m')

  if args.use_pdb:
    new_pdb_ids = pd.read_csv(args.use_pdb, comment='#', header=None, sep='\s+').iloc[:,:].values.tolist()
    print('\n  \033[34m## Work with a list of pre-supplied PDB: \033[31m{}\033[0m'.format(len(new_pdb_ids)))
  else:
    new_pdb_ids = None

#################
  ## Get basic settings from parameter file
  work_dir   = parms['WORKDIR']
  outpref    = parms['OUTPREF']
  len_cutoff = int(parms['LENCUTOFF'])
  idt_cutoff = float(parms['IDTCUTOFF'])
  known_pdbs = ReadKnownKinaseList(parms['RDKNOWNKIN'])
  read_nonkinase = parms['RDNONKIN']
  read_checked   = parms['RDCHECKED']

################

  ## Read from a list of PDBs for known non-kinase, will be used to 
  ## subtract from retrieved entries
  if not re.search('none', read_nonkinase, re.IGNORECASE):
    if os.path.isfile(read_nonkinase):
      non_kinase = pd.read_csv(read_nonkinase,sep='\s+',comment='#',header=None,names=['pdb_idx'])
    else:
      print('  \033[31m#1# Param ERROR: Supplied Non-kinase List file is not found: \033[0m\n'+read_nonkinase)
      non_kinase = pd.DataFrame([], columns=['pdb_idx'])
  else:
    non_kinase = pd.DataFrame([], columns=['pdb_idx'])
  print('  > Found Known Non-Kinase entries: \033[31m{}\033[0m'.format(len(non_kinase)))


  ## Look for previously checked PDBs from the query list
  ## to be used to Subtract reviewed PDB from the query list
  if not re.search('none', read_checked, re.IGNORECASE):
    if os.path.isfile(read_checked):
      checked_df = pd.read_csv(read_checked,sep='\s+',header=None,comment='#',names=['pdb_id'])
    else:
      print('  \033[31m#1# Param ERROR: Supplied Checked-Kinase List file is not found: \033[0m\n'+read_checked)
      checked_df = pd.DataFrame([], columns=['pdb_id'])
  else:
    checked_df = pd.DataFrame([], columns=['pdb_id'])
  print('\n  > Found checked-kinase entries: \033[31m{}\033[0m'.format(len(checked_df)))


################################################################
  new_ents = []  # store new pdb entries found from queries

  ## Skip all the checking steps and directly download/extract new PDB if 'use_pdb'
  if new_pdb_ids is None:
    ## create the predefined queries to search for all protein kinases
    kinase_queries = KinaseDownloadQueries( parms )

    ## Get the list of PDBs that hits the queries, filter out the known ones
    print('\n  \033[34m## Checking RCSB PDB REST API for entries hitting the "Kinase" queries ##\033[0m')
    query_rt = [SearchRCSBWithQuery(key, val) for key, val in kinase_queries.items()]
    query_et = list(set([pdb_id for sublist in query_rt for pdb_id in sublist]))  # unpack
    que_ents = sorted(CheckExistingPDBs( query_et, known_pdbs ))
    print('\n  > Found new PDB entries: \033[31m{}\033[0m'.format(len(que_ents)))

    ## Subtract previously checked PDB (kinase or nonkinase) from the entry list
    if not re.search('none', read_checked, re.IGNORECASE):
      print('  > Found Previously-Checked-PDB List: \033[31m{}\033[0m'.format(read_checked))
      que_df   = pd.DataFrame(que_ents, columns=['pdb_id'])
      new_ents = que_df[~que_df.pdb_id.isin(checked_df.pdb_id.tolist())].pdb_id.tolist()
      print('  > Entries after subtracting from known non-kinases: \033[31m{}\033[0m'.format(len(new_ents)))
    else:
      new_ents = que_ents

    ## End the job if no new PDB is found
    if not len(new_ents):
      sys.exit('\n  \033[31m### Checked and No New PDB Found ###\033[0m\n')


#####################################
    ## download FASTA of the new entries
    print('\n  \033[34m## Downloading FASTA of the New Entries ##\033[0m')
    dl_fasta = [ DownloadFASTA(pdb_id) for pdb_id in tqdm(new_ents) ]
    raw_fasta = list([pdb_id for sublist in dl_fasta for pdb_id in sublist])  # unpack
    print('  > FASTA search of entries: {0} -> \033[31m{1}\033[0m chains'.format(len(new_ents), len(raw_fasta)))
    raw_f_df = pd.DataFrame(raw_fasta, columns=['pdb_idx','pdb_id','chain_id','length','seq'])

    ## Subtract retrieved entries that are known to be non-kinase
    if not re.search('none', read_nonkinase, re.IGNORECASE):
      raw_f_df = raw_f_df[~raw_f_df.pdb_idx.isin(non_kinase.pdb_idx.tolist())].reset_index()
      print('  > Entries after subtracting from known non-kinases: \033[31m{}\033[0m'.format(len(raw_f_df)))

    ## remove sequences too short to be kinase catalytic domain
    ## remove duplicate sequences of different chains from same pdb_id, keep first one
    clean_f_df = raw_f_df[ raw_f_df.length > len_cutoff ]
    uniq_f_df  = clean_f_df.drop_duplicates(subset=['pdb_id','seq'], keep='first').reset_index(drop=True)

    print('\n  \033[34m## Removing FASTA \033[31m < {}\033[0m \033[34mresidues ##\033[0m'.format(len_cutoff))
    print('   > \033[31m{0}\033[0m/{1} FASTA passed the filter'.format(len(clean_f_df),len(raw_f_df)))
    print('   > \033[31m{0}\033[0m/{1} FASTA are unique'.format(len(uniq_f_df),len(clean_f_df)))


#####################################
    ## Check if FASTA is kinase by comparing sequence identity to all canonical human kinases (kinome)
    ## output is Length>175, SeqIdent>40%, [[pdb_id, chain_id], [], [] ]
    new_pdb_ids, nonkin_seq_df = CheckKinaseSeqIdentity( uniq_f_df, parms['KINOMEDB'], len_cutoff, idt_cutoff, outpref )
    print('  > Found kinase seq: \033[31m{0}/{1}\033[0m'.format(len(new_pdb_ids), len(uniq_f_df)))

    ## collect all non_kinases and write it out for future reference
    nonkin_df = pd.concat( [non_kinase, raw_f_df[ raw_f_df.length <= len_cutoff ], nonkin_seq_df ], sort=True).drop_duplicates(subset='pdb_idx')
    nonkin_df.pdb_id = nonkin_df.pdb_idx.apply(lambda x: x.split('_')[0])
    nonkin_df.to_csv(parms['WTNONKIN'], header=None, index=None, columns=['pdb_idx'],sep='\t')
    print('  > Found non-kinase PDB: \033[31m{}\033[0m'.format(len(nonkin_df)))
    print('  > Saved non-kinase PDB: '+parms['WTNONKIN'])


########### new_pdb_ids #############
  ## By this stage, if 'new_pdb_ids' is still None, meaning no new PDB is found. Exit
  if new_pdb_ids is None:
    sys.exit('\n  \033[31m### Checked and No New PDB Found ###\033[0m\n')


#####################################
  mpi = multiprocessing.Pool()
  ## Download new PDB, extract the chain(s)
  print('\n  \033[34m## Downloading new PDB ##\033[0m')
  oPD  = DownloadNewPDB(work_dir=work_dir)
  Temp = [x for x in tqdm(mpi.imap(oPD, new_pdb_ids), total=len(new_pdb_ids))]
#  Temp = [ oPD(info) for info in tqdm(new_pdb_ids) ]
  raw_pdbs = [pdb for pdb in Temp if pdb is not None]
  print('  > Downloaded new PDB entries: \033[31m{0}/{1}\033[0m'.format(len(raw_pdbs),len(new_pdb_ids)))

  ## updated known kinase list by adding new ones to it; chain_id is also recorded
  for item in new_pdb_ids:
    if item[0] not in known_pdbs:
      known_pdbs[item[0]] = [item[1]]
    else:
      known_pdbs[item[0]].append(item[1])

  ## write out all updated + exisitng kinase structures (with output extension)
  print('\n  \033[34m## Finalize all results ##\033[0m')
  with open(parms['WTNEWKNOWN'], 'w') as fo:
    for i in sorted(known_pdbs.keys()):
      for chain in known_pdbs[i]:
        fo.write('{0}_{1}.{2}\n'.format(i, chain, parms['OUTEXT']))
  print('\n  > Total number of known kinase PDB_x structures: \033[31m{}\033[0m'.format(len(known_pdbs)))

  ## write out all checked PDB, may use this next time to skip these checked ones
  with open(parms['WTCHECKED'], 'w') as fc:
    for pdb_id in sorted(list(set(checked_df.pdb_id.tolist() + new_ents + list(known_pdbs.keys() ))) ):
      fc.write(pdb_id+'\n')
  print('\n  \033[34m> All checked PDBs are listed in:\033[0m \033[31m{}\033[0m'.format(parms['WTCHECKED']))


#####################################
  ## Generate MSA alignment for new PDBs against known MSA alignment
  refname   = parms['REFPDB'].split('/')[-1].split('.pdb')[0]

  print('\n  \033[34m## Get MSA alignment for PDB structures ##\033[0m')
  
  # Use canonical kinome dbase to align; suppress Blastp StdOut during MSA alignment
  Fasta_Dict = CacheSeqDatabase(parms['KINOMEGPDB'])
  oAl = GenerateProfileAlignment( hom_dir=parms['WORKDIR'],  tmp_dir=parms['WORKDIR'], 
                                  rst_dir=parms['WORKDIR'],  ref_pdb=refname, 
                                  f_nogap=parms['KINOMEDB'], f_dict=Fasta_Dict )    

  Aligned_Fasta = [x for x in tqdm(mpi.imap(oAl, raw_pdbs),total=len(raw_pdbs))]
#  Aligned_Fasta = [ oAl(pdb) for pdb in tqdm(raw_pdbs) ]
  mpi.close()
  mpi.join()

  # Remove any input kinase that have been rejected due to low seq identity
  with open('{0}.fasta'.format(outpref), 'w') as fo:
    SeqIO.write( [Fasta_Dict[refname]]+[x for x in Aligned_Fasta if x is not None], fo, 'fasta' )


#####################################
  ## Superpose new PDBs to 1ATP C-lobe
  print('\n  \033[34m## Superposing all PDB structure with extension \033[31m{}\033[0m \033[34m##\033[0m'.format(parms['OUTEXT']))
  SuperposePDB( parms['PYMOL'],  parms['REFPDB'], raw_pdbs, None, 
                parms['OUTEXT'], parms['REFRES'], outpref )


##########################################################################################
def cmdlineparse():
  p = ArgumentParser(description="command line arguments")

  p.add_argument('-set', action='store_true', help='Print out Parameter file used by "-r"')
  p.add_argument('-r', dest='read_param', required=True,
                  help='Read in Parameter file')
  p.add_argument('-p', dest='use_pdb', required=False,
                  help='List of PDB (and Chain ID) to directly download (format: pdb_id\tchain_id)')

  return p.parse_args()

########################################################################
if __name__ == '__main__':
  main()

#
#  v1.0  20.01.10
#
#