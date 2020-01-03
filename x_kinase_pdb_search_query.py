#!/usr/bin/env python3

import os
import re
import sys
import time
import requests
import pandas as pd

from Bio import PDB

############################################################################################
## search RCSB PDB site with predefined kinase search queries, 
## download those that are not in internal library already
def DownloadFASTA( pdb_id ):
  peptides = []

  ## use 'requests' to download fasta and output the fasta
  download_fasta = 'https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList={0}&compressionType=uncompressed'.format(pdb_id)
  fasta_url = requests.get(download_fasta)

  ## delay retrieve by 0.1s to avoid server from download lockout (prevent downloading)
  time.sleep(0.1)

  ## downloaded FASTA is byte, need to split it into individual chains by '>',
  ## then concatenate pdb_id and chain_id. Split string by 'SEQUENCE', then
  ## remove all '\n'
  Tmp = [ x for x in fasta_url.content.decode().split('>') if x is not '' ]
  for line in Tmp:
    try:
      pdb_tmp, seq_tmp = re.sub(':', '_', line).split('SEQUENCE')
    except ValueError:
      print(' > \033[31mFailed to retrieve:\033[0m '+pdb_id)
      continue
    pdb_idx = pdb_tmp.split('|')[0]
    pdb_id, chain_id = pdb_idx.split('_')
    seq     = re.sub('\n', '', seq_tmp)
    length  = len(seq)
    peptides.append([pdb_idx, pdb_id, chain_id, length ,seq])

  return peptides


############################################################################################
## Get all PDB IDs that hit the query, compare query PDBs to known list and
## output the new ones not in the known list
def SearchRCSBWithQuery( key, search_query ):
  pdb_url = 'http://www.rcsb.org/pdb/rest/search'
  header = {'Content-Type': 'application/x-www-form-urlencoded'}

  ## seach RCSB with query using requests
  response = requests.post(pdb_url, data=search_query, headers=header)
  if response.status_code == 200:
    print("  > Found PDB entries matching query: \033[31m{} - {}\033[0m".format(key, len(response.text)))
  else:
    sys.exit("\033[34m FATAL \033[0m- Failed to retrieve results")

  ## response.content will give a list of 4-letter PDB_ID
  pdb_id_list = sorted([i for i in response.content.decode().split('\n') if i != ''])
  print('  > Parsed PDB entries that matches: \033[31m{}\033[0m'.format(len(pdb_id_list)))

  return pdb_id_list


#############################################################################
## Check if the newly retreived PDB name already exists in internal library
def ReadKnownKinaseList( list_name ):
  known_pdbs = {}
  with open(list_name, 'r') as fi:
    for l in fi:
      if l.rstrip() and not re.search('#', l):
        pdb_id, rest = l.rstrip().split('_')
        chain_id     = rest.split('.')[0]
        if pdb_id not in known_pdbs:
          known_pdbs[pdb_id] = [chain_id]
        else:
          known_pdbs[pdb_id].append(chain_id)
  return known_pdbs

## First open the known_pdb.list, check if query PDB is in known_pdb.list already
## return those that are not found in known_pdb.list
def CheckExistingPDBs( pdb_list, known_pdbs ):
  return [ pdb_id for pdb_id in pdb_list if pdb_id not in known_pdbs ]


##########################################################################################
## Download PDB from RCSB then extract the Chain of interest
def DownloadNewPDB( work_dir, info ):
  pdb_id, chain_id = info

  ## BioPython downloads PDB but it gives a lowercase name in pdb{}.ent format
  PDB.PDBList(verbose=False).retrieve_pdb_file(pdb_id, pdir=work_dir, obsolete=False, file_format='pdb')
  biopdb_name = '{0}/pdb{1}.ent'.format(work_dir, pdb_id.lower())

  ## Read the PDB file and extract the chain from structure[0]
  model = PDB.PDBParser(PERMISSIVE=1,QUIET=True).get_structure(pdb_id, biopdb_name)[0]

  io = PDB.PDBIO()
  io.set_structure(model[chain_id])
  io.save('{0}/{1}_{2}.pdb'.format(work_dir, pdb_id, chain_id))
  os.system('mv {1} {0}/{2}.ent ; bzip2 -f {0}/{2}.ent'.format(work_dir, biopdb_name, pdb_id))

  return '{0}/{1}_{2}.pdb'.format(work_dir, pdb_id, chain_id)


############################################################################################
  ## Check FASTA is kinase by comparing sequence identity to all human kinases
  ## output is a Dictionary of pdb_id to chain_id
def CheckKinaseSeqIdentity( uniq_f_df, kinase_db, len_cutoff, idt_cutoff, outpref ):

  # good:  length match > 175, seq ident > 40%
  # check: length match > 175, seq ident 30-40%
  # bad:   length match > 175, seq ident < 30%
  # no:    length match <= 175
  f_ok  = open(outpref+'.good_seq_ident.txt', 'w')
  f_ok.write('pdb_idx\tlength\tident\tsimi\n')
  f_chk = open(outpref+'.check_seq_ident.txt', 'w')
  f_chk.write('pdb_idx\tlength\tident\tsimi\n')
  f_bad = open(outpref+'.bad_seq_ident.txt', 'w')
  f_bad.write('pdb_idx\tlength\tident\tsimi\n')
  f_no  = open(outpref+'.no_seq_ident.txt', 'w')
  f_no.write('pdb_idx\tlength\tident\tsimi\n')

  new_pdb_ids = []    # collection of confirmed kinases
  non_kin_idx = []    # collection of non-kinase pdb_idx
  print('\n \033[34m## Comparing sequence identity of unique PDB FASTA ##\033[0m]')
  for idx, x in tqdm(uniq_f_df.iterrows(), total=len(uniq_f_df)):
    fasta_file = '{0}_{1}.fasta'.format(x.pdb_id, x.chain_id)
    with open(fasta_file, 'w') as fo:
      fo.write('>'+x.pdb_idx+'\n'+x.seq+'\n')

    Data = BlastpPairwiseIdentity( '.', fasta_file, kinase_db )
    if Data is None:
      f_no.write(fasta_file+'\t no output\n')
      continue
    imat_df = pd.DataFrame(Data, columns=['kinase','length','ident','simi'])

    ## Hard cutoff of required matching to 175 residues, if fewer than 200 res,
    ## unlikely to be a kinase catalytic domain. if seq ident < 40% of any 
    ## known human kinases , need to check and confirm if it is a bacterial
    ## or viral kinases; mammalian kinases are very similar, even chicken/fish
    sele_df = imat_df[ imat_df.length > (len_cutoff - 40) ]
    info = '{0}_{1}\t{2}\t{3}\t{4}\n'.format(x.pdb_id, x.chain_id, 
                imat_df.length[0], imat_df.ident[0], imat_df.simi[0])

    ## catch a strange error with ident.iloc[0] key error
    try:
      test = sele_df.ident
    except KeyError:
      print('  \033[31mERROR:\033[0m "sele_df.ident[0]" - '+x.pdb_idx)
      continue

    if len(sele_df) == 0 or len(sele_df.ident) == 0:
      non_kin_idx.append('{0}_{1}'.format(x.pdb_id, x.chain_id))
      f_no.write(info)
      continue
    elif sele_df.ident.iloc[0] > idt_cutoff:
      f_ok.write(info)  
      new_pdb_ids.append( [x.pdb_id, x.chain_id] )
    elif sele_df.ident.iloc[0] > (idt_cutoff - 10.):
      f_chk.write(info)
    else:
      non_kin_idx.append('{0}_{1}'.format(x.pdb_id, x.chain_id))
      f_bad.write(info)

  f_bad.close()
  f_chk.close()
  f_ok.close()
  f_no.close()

  nonkin_seq_df = uniq_f_df[ uniq_f_df.pdb_idx.isin(non_kin_idx) ]

  return new_pdb_ids, nonkin_seq_df


########################################################################


########### These Queries are generalized Keyword search #############
## Search RCSB Protein Data Bank using the 'AdvancedKeywordQuery' keywords:, then
## select for 'X-RAY structures' -- NMR are unlikely to be the catalytic domain
## These queries are derived from the 'Query Details' but replacing some
## of the key factors


### EnzymeClassificationTree Search for 2.7.11.1: Non-specific serine/threonine protein kinase
### Experimental Method is X-RAY
n_stk_query = '''
<orgPdbCompositeQuery version="1.0">
 <queryRefinement>
  <queryRefinementLevel>0</queryRefinementLevel>
  <orgPdbQuery>
    <queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>
    <description>Text Search for: 2.7.11.1</description>
    <keywords>2.7.11.1</keywords>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>1</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <queryType>org.pdb.query.simple.ExpTypeQuery</queryType>
    <description>Experimental Method is X-RAY</description>
    <mvStructure.expMethod.value>X-RAY</mvStructure.expMethod.value>
    <mvStructure.expMethod.exclusive>y</mvStructure.expMethod.exclusive>
  </orgPdbQuery>
 </queryRefinement>
</orgPdbCompositeQuery>
'''

### EnzymeClassificationTree Search for 2.7.11.30: Receptor protein serine/threonine kinase	
### Experimental Method is X-RAY	
r_stk_query = '''
<orgPdbCompositeQuery version="1.0">
 <queryRefinement>
  <queryRefinementLevel>0</queryRefinementLevel>
  <orgPdbQuery>
    <queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>
    <description>Text Search for: 2.7.11.30</description>
    <keywords>2.7.11.30</keywords>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>1</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <queryType>org.pdb.query.simple.ExpTypeQuery</queryType>
    <description>Experimental Method is X-RAY</description>
    <mvStructure.expMethod.value>X-RAY</mvStructure.expMethod.value>
    <mvStructure.expMethod.exclusive>y</mvStructure.expMethod.exclusive>
  </orgPdbQuery>
 </queryRefinement>
</orgPdbCompositeQuery>
'''

### EnzymeClassificationTree Search for 2.7.10.2: Non-specific protein-tyrosine kinase
### Experimental Method is X-RAY	
n_tk_query = '''
<orgPdbCompositeQuery version="1.0">
 <queryRefinement>
  <queryRefinementLevel>0</queryRefinementLevel>
  <orgPdbQuery>
    <queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>
    <description>Text Search for: 2.7.10.2</description>
    <keywords>2.7.10.2</keywords>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>1</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <queryType>org.pdb.query.simple.ExpTypeQuery</queryType>
    <description>Experimental Method is X-RAY</description>
    <mvStructure.expMethod.value>X-RAY</mvStructure.expMethod.value>
    <mvStructure.expMethod.exclusive>y</mvStructure.expMethod.exclusive>
  </orgPdbQuery>
 </queryRefinement>
</orgPdbCompositeQuery>
'''

### EnzymeClassificationTree Search for 2.7.10.1: Receptor protein-tyrosine kinase
### Experimental Method is X-RAY
r_tk_query = '''
<orgPdbCompositeQuery version="1.0">
 <queryRefinement>
  <queryRefinementLevel>0</queryRefinementLevel>
  <orgPdbQuery>
    <queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>
    <description>Text Search for: 2.7.10.1</description>
    <keywords>2.7.10.1</keywords>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>1</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <queryType>org.pdb.query.simple.ExpTypeQuery</queryType>
    <description>Experimental Method is X-RAY</description>
    <mvStructure.expMethod.value>X-RAY</mvStructure.expMethod.value>
    <mvStructure.expMethod.exclusive>y</mvStructure.expMethod.exclusive>
  </orgPdbQuery>
 </queryRefinement>
</orgPdbCompositeQuery>
'''

### EnzymeClassificationTree Search for 2.7.12.1: Dual-specificity kinase	
### Experimental Method is X-RAY	
dual_k_query = '''
<orgPdbCompositeQuery version="1.0">
 <queryRefinement>
  <queryRefinementLevel>0</queryRefinementLevel>
  <orgPdbQuery>
    <queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>
    <description>Text Search for: 2.7.12.1</description>
    <keywords>2.7.12.1</keywords>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>1</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <queryType>org.pdb.query.simple.ExpTypeQuery</queryType>
    <description>Experimental Method is X-RAY</description>
    <mvStructure.expMethod.value>X-RAY</mvStructure.expMethod.value>
    <mvStructure.expMethod.exclusive>y</mvStructure.expMethod.exclusive>
  </orgPdbQuery>
 </queryRefinement>
</orgPdbCompositeQuery>
'''

def KinaseDownloadQueries( parms ):

  return {
    'NonSpec_stk': n_stk_query,    # non-specific S/T kinase
    'Recp_stk':    r_stk_query,    # receptor S/T kinase
    'NonSpec_tk':  n_tk_query,     # non-specific Tyr kinase
    'Recp_tk':     r_tk_query,     # receptor Tyr kinase
    'DualSpec_k':  dual_k_query    # dual-specific kinase
    }