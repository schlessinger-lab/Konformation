#!/usr/bin/env python3

import io
import os
import re
import sys
import time
import requests
import pandas as pd

from tqdm import tqdm
from pathos import multiprocessing

from Bio import PDB
from Bio.PDB import Select

from x_search_align import BlastpPairwiseIdentity
from x_pdb_modif_resid import ReplacePDBModifiedAA

#######################################################################################
#
#  v1.0  20.01.02
#
# Bug - 20.01.12 - Biopython 1.74/1.76 have issue with PDB writeout. v1.72 has no issue
#   > AttributeError: 'Atom' object has no attribute 'disordered_get_list'
#   to go around this issue, generate an intermediate file instead
#
#
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
    print("  > Found PDB entries matching query: \033[31m{0} - {1}\033[0m".format(key, len(response.text)))
  else:
    sys.exit("\033[34m FATAL \033[0m- Failed to retrieve results")

  ## response.content will give a list of 4-letter PDB_ID
  pdb_id_list = sorted([i for i in response.content.decode().split('\n') if i != ''])
  print('  > Parsed PDB entries that matches: \033[31m{0}\033[0m'.format(len(pdb_id_list)))

  return pdb_id_list


#############################################################################
## Check if the newly retreived PDB name already exists in internal library
def ReadKnownKinaseList( list_name ):
  known_pdbs = {}
  if list_name == 'None':
    return known_pdbs

  with open(list_name, 'r') as fi:
    for l in fi:
      if l.rstrip() and not re.search('#', l):
        pdb_id, rest = l.rstrip().split('_')
        chain_id     = rest.split('.')[0]
        if pdb_id not in known_pdbs:
          known_pdbs[pdb_id] = [chain_id]
        else:
          known_pdbs[pdb_id].append(chain_id)
  print('  > Found Known Kinase entries: \033[31m{0}\n - {1}\033[0m\n'.format(len(known_pdbs), list_name))
  return known_pdbs

## First open the known_pdb.list, check if query PDB is in known_pdb.list already
## return those that are not found in known_pdb.list
def CheckExistingPDBs( pdb_list, known_pdbs ):
  return [ pdb_id for pdb_id in pdb_list if pdb_id not in known_pdbs ]


##########################################################################################
## Download PDB from RCSB then extract the Chain of interest, return path/name of PDB
class DownloadNewPDB(object):
  def __init__( self, work_dir='' ):
    self.work_dir = work_dir
  def __call__( self, info ):
    return self.download_pdb( info )

  def download_pdb( self, info ):
    pdb_id, chain_id = info

    ## Check if atom has alternative position, if so, keep 'A' position and remove the flag
    ## but somehow this class doesn't seem to function well
    class NotDisordered(Select):
      def accept_atom(self, atom):
        if not atom.is_disordered() or atom.get_altloc() == 'A':
          atom.set_altloc(' ')
          return True
        else:
          return False

    ## BioPython downloads PDB but it gives a lowercase name in pdb{}.ent format
    biopdb_name = '{0}/pdb{1}.ent'.format(self.work_dir, pdb_id.lower())
    biopdb_modf = '{0}/pdb{1}.mod.ent'.format(self.work_dir, pdb_id.lower())
    if not os.path.isfile(biopdb_modf):
      try:
        PDB.PDBList(verbose=False).retrieve_pdb_file(pdb_id, pdir=self.work_dir, 
                                                obsolete=False, file_format='pdb')
      except FileNotFoundError:
        print('  \033[31m> ERROR: BioPython cannot download PDB: \033[0m'+pdb_id)
        return None

    ## Replace modified AA to avoid mis-recognition in biopython readin
    ## Replace disordered atoms and keep only the "A" variant
    ReplacePDBModifiedAA(biopdb_name, biopdb_modf)
    os.system('grep "REMARK  " {0} > {0}.remark'.format(biopdb_modf))
    with open(biopdb_modf, 'r') as fi:
      remarks = [l for l in fi if re.search('REMARK HET ',l)]

    ## Read the PDB file and extract the chain from structure[0]
    try:
      model = PDB.PDBParser(PERMISSIVE=1,QUIET=1).get_structure(pdb_id, biopdb_modf)[0]
    except KeyError:
      print('  \033[31m> ERROR: BioPython cannot read in PDB: \033[0m'+biopdb_modf)
      return None
    except ValueError:
      print('  \033[31m> ERROR: PDB file is empty: \033[0m'+biopdb_modf)
      return None

    ### Bug alert: as of 20.02.18, Biopython dev hasn't come up with good
    ### strategy to fix the 'atom.disordered_get_list()' issue with alternative
    ### position of residue side chains. To go around this, will physically
    ### remove "B" variant and keep only "A" variant in 
    io = PDB.PDBIO()
    io.set_structure(model[chain_id])
    io.save('{0}/{1}_{2}.pdb'.format(self.work_dir, pdb_id, chain_id), 
            select=NotDisordered())

    # Attach REMARK to end of PDB as safekeeping
    os.system('cat {0}/{1}_{2}.pdb {3}.remark > {1}.temp'.format(self.work_dir, pdb_id, chain_id, biopdb_modf))
    os.system('mv {1}.temp {0}/{1}_{2}.pdb'.format(self.work_dir, pdb_id, chain_id))
#    os.system('mv {1} {0}/{2}.ent'.format(self.work_dir, biopdb_name, pdb_id))
#    os.system('bzip2 -f {0}/{1}.ent'.format(self.work_dir, pdb_id))
#    os.system('rm {0} {0}.remark'.format(biopdb_modf))

    return '{0}/{1}_{2}.pdb'.format(self.work_dir, pdb_id, chain_id)


############################################################################################
def RunBlastDB( data ):
  kinase_db, pdb_id, chain_id, pdb_idx, seq = data

  fasta_file = '{0}_{1}.fasta'.format(pdb_id, chain_id)
  with open(fasta_file, 'w') as fo:
    fo.write('>'+pdb_idx+'\n'+seq+'\n')

  Identity = BlastpPairwiseIdentity( '.', fasta_file, kinase_db )
  if Identity is not None:
    imat_df = pd.DataFrame(Identity, columns=['kinase','length','ident','simi'])
  else:
    imat_df = None

  return [data, imat_df]


############################################################################################
## Check FASTA is kinase by comparing sequence identity to all canonical human kinases (kinome)
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
  print('\n \033[34m## Comparing sequence identity of unique PDB FASTA ##\033[0m')

  ## Run Blastp in parallel
  Data = [[kinase_db, r.pdb_id, r.chain_id, r.pdb_idx, r.seq] for idx, r in uniq_f_df.iterrows()]
  mpi  = multiprocessing.Pool()
  df_l = [x for x in tqdm(mpi.imap(RunBlastDB, Data), total=len(Data))]
  mpi.close()
  mpi.join()

  for items in df_l:
    data, imat_df = items
    kinase_db, pdb_id, chain_id, pdb_idx, seq = data
    fasta_file = '{0}_{1}.fasta'.format(pdb_id, chain_id)

    if imat_df is None:
      f_no.write(fasta_file+'\t no output\n')
      continue

    ## Required matching to at least 200 residues as cutoff, if fewer than 200 res,
    ## unlikely to be a kinase catalytic domain. if seq ident < 40% of any 
    ## known human kinases, need to check and confirm if it is a bacterial
    ## or viral kinases; mammalian kinases are very similar, even to chicken/fish
    sele_df = imat_df[ imat_df.length > (len_cutoff - 40) ]
    info = '{0}_{1}\t{2}\t{3}\t{4}\n'.format(pdb_id, chain_id, 
                imat_df.length[0], imat_df.ident[0], imat_df.simi[0])

    ## catch a strange error with ident.iloc[0] key error
    try:
      test = sele_df.ident
    except KeyError:
      print('  \033[31mERROR:\033[0m "sele_df.ident[0]" - '+pdb_idx)
      continue

    if len(sele_df) == 0 or len(sele_df.ident) == 0:
      non_kin_idx.append('{0}_{1}'.format(pdb_id, chain_id))
      f_no.write(info)
      continue
    elif sele_df.ident.iloc[0] > idt_cutoff:
      f_ok.write(info)  
      new_pdb_ids.append( [pdb_id, chain_id] )
    elif sele_df.ident.iloc[0] > (idt_cutoff - 10.):
      f_chk.write(info)
    else:
      non_kin_idx.append('{0}_{1}'.format(pdb_id, chain_id))
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