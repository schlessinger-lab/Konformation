#!/usr/bin/env python3

import re
import os
import sys
import time
import requests
import xmltodict

import numpy as np
import pandas as pd

from aa_residue import AA
from aa_residue import UnnaturalAA
from aa_residue import SaltAdditive

from tqdm import tqdm

##########################################################################   
## Read in pdb list with different formats:
# Format of the list:
#  1) <pdb>_<chain_id>.xxx.pdb
#  2) <pdb>_<chain_id>.xxx.pdb <chain_id>
#  3) <pdb>.xxx.pdb            <chain_id>
#  4) <pdb>                    <chain_id>
#  (PDB ID and chain ID separated by '_')
#
#  e.g.:   1ATP_E.xxx.pdb
#          1ATP_E
#          3HHP.xxx.pdb    C
#          6GTT            A
#          6GTT_A          A  cidi
#          3HHP_A          C
def ProcessPDBList( pdb_list ):
  PDB_List = pd.read_csv(pdb_list, sep='\s+', header=None, comment='#').to_numpy()

  PDB   = []
  for entry in PDB_List:
    conf = None
    if len(entry) == 1:
      if re.search(r'.pdb', entry[0]):
        pdb_id   = entry[0].split('.')[0]
        pdb      = entry[0].split('.')[0].split('_')[0]
        chain_id = entry[0].split('.')[0].split('_')[1]
      else:
        pdb_id   = entry[0]
        pdb      = entry[0].split('_')[0]
        chain_id = entry[0].split('_')[1]
      PDB.append([pdb_id, pdb, chain_id, conf])

    elif len(entry) >= 2:
      if re.search(r'.pdb', entry[0]):
        if re.search('_', entry[0]):
          pdb_id   = entry[0].split('.')[0]
          pdb      = entry[0].split('.')[0].split('_')[0]
          chain_id = entry[0].split('.')[0].split('_')[1]
        else:
          pdb      = entry[0].split('.')[0]
          chain_id = entry[1]
          pdb_id   = pdb+'_'+chain_id
      else:
        if re.search('_', entry[0]):
          pdb_id   = entry[0]
          pdb      = entry[0].split('_')[0]
          chain_id = entry[0].split('_')[1]
        else:
          pdb      = entry[0]
          chain_id = entry[1]
          pdb_id   = pdb+'_'+chain_id
      if len(entry) == 3:
        conf = entry[2]
      PDB.append([pdb_id, pdb, chain_id, conf])
    else:
      sys.exit('\033[31m ERROR: PDB list format is invalid\033[0m')

  return PDB


##########################################################################
#   Format for -c <csv> header:
#   'pdb_id','Class','cidi_prob','cido_prob','codi_prob',...
def ProcessCSVList( csv_list ):

  def get_pdb_id( inp ):
    return inp.split('_')[0]
  def get_chain_id( inp ):
    return inp.split('_')[1]

  PDB_df = pd.read_csv(csv_list, sep=',', comment='#')
  PDB_List = PDB_df[['pdb_id','Class']]

  PDB_List['pdb']      = PDB_List['pdb_id'].apply(get_pdb_id)
  PDB_List['chain_id'] = PDB_List['pdb_id'].apply(get_chain_id)

  return PDB_List[['pdb_id','pdb','chain_id','Class']].to_numpy()


##########################################################################
## Download the PDB page from the Web and search for the UniProt ID
def SearchPDB( inp ):
  pdb_id, pdb, chain_id, conf = inp

  pdb_length, uni_id, tax_id, species = None, None, None, None
  ec, pmid, p_name, uni_id = None, None, None, None
  resolu, deposit, release, latest = None, None, None, None

  ## Downloaded PDB info, convert XML into dict and extract data
  html_pdb = 'https://www.rcsb.org/pdb/rest/describePDB?structureId={0}'.format(pdb)
  pdb_info = requests.get(html_pdb)
  pdb_dict = xmltodict.parse(pdb_info.content.decode())
  pdbx     = pdb_dict['PDBdescription']['PDB']

  ## not all values are available, some can be missing
  if '@pubmedId' in pdbx:         pmid   = pdbx['@pubmedId']
  if '@resolution' in pdbx:       resolu = pdbx['@resolution']
  if '@deposition_date' in pdbx:  deposit= pdbx['@deposition_date']
  if '@release_date' in pdbx:     release= pdbx['@release_date']
  if '@last_modification_date' in pdbx: latest = pdbx['@last_modification_date']
  if '@citation_authors' in pdbx: authors= pdbx['@citation_authors']

## delay retrieve by 0.1s to avoid server from download lockout (prevent downloading)
  time.sleep(0.1)

  ## Download individual PDB chain info, convert XML into dict and extract
  html_chain = 'https://www.rcsb.org/pdb/rest/describeMol?structureId={0}.{1}'.format(pdb, chain_id)
  chain_info = requests.get(html_chain)
  chain_dict = xmltodict.parse(chain_info.content.decode())
  chnx       = chain_dict['molDescription']['structureId']['polymer']

  if '@length' in chnx:   pdb_length = chnx['@length']
  if '@weight' in chnx:   weight     = float(chnx['@weight'])
  if 'Taxonomy' in chnx:
    if '@id' in chnx['Taxonomy']:   tax_id     = chnx['Taxonomy']['@id']
    if '@name' in chnx['Taxonomy']: species    = chnx['Taxonomy']['@name']
  if 'macroMolecule' in chnx:
    if '@name' in chnx['macroMolecule']: p_name     = chnx['macroMolecule']['@name']
    if 'accession' in chnx['macroMolecule']:
      if '@id' in chnx['macroMolecule']['accession']:
        uni_id     = chnx['macroMolecule']['accession']['@id']
  if 'enzClass' in chnx:  ec         = chnx['enzClass']['@ec']


#########################################
  ## Find UniProt ID, chain ID, chain Length, in PDB webpage

  Entity = {}
  Entity[chain_id] = HeadData(pdb_id=pdb_id, pdb=pdb, chain_id=chain_id, 
                              conf=conf, pdb_length=pdb_length, uni_id=uni_id,
                              ec=ec, taxid=tax_id, species=species, 
                              p_name=p_name, pmid=pmid, resolu=resolu, 
                              latest=latest, deposit=deposit, release=release)
  Entity[chain_id].chain = chain_id

################################
  # Retrieve data from HEADER data in RCSB webpages. Data is not XML format
  html_head = 'https://files.rcsb.org/header/{0}.pdb'.format(pdb)
  head_info = requests.get(html_head)
  headx = [ x for x in head_info.content.decode().split('\n') if x is not '' ]

  HD = {}
  mol_id, xtalc_remark = 0, 0
  space, mutation, aa_modif, ligand, salt = None, None, None, None, None
  XtalC, Het_Lns, HetNam_Lns, ModRes, HetNam = [], [], [], [], []
  SeqAdv, DBRef = [], []

  ## Find Resolution, publish date, packing group in PDB header
  ## These data only needed once and independent of chain_id, put them into 
  ## the object of each individual chain_id found in PDB header
  for l in headx:
    # Get crystal unit cell
    if re.search(r'CRYST', l):
      crystal = l.rstrip()
    # Get Mutation and protein-uniprot data, if available
    # ideally get them to sync up to chain ID but too hard to parse
    if re.search(r'^SEQADV', l):
      SeqAdv.append(l.rstrip())
    if re.search(r'^DBREF', l):
      DBRef.append(l.rstrip())

    if re.search(r'^REMARK', l):
      # Get space group
      if re.search(r'SPACE GROUP:', l):
        space = l.rstrip().split(':')[1].strip()
      # Get crystalization condition
      if re.search(r'CRYSTALLIZATION', l):
        XtalC.append(l)
        xtalc_remark = l.split()[1]
      if l.split()[1] == xtalc_remark:
        XtalC.append(l)

    # Get modified residues
    if re.search(r'^MODRES', l):
      ModRes.append(l.rstrip().split()[2])
#      print('MODRES: '+l.rstrip().split()[2])
    # Get ligand PDB_ID
    if re.search(r'^HET ', l):
      Het_Lns.append(l.rstrip().split()[1])
    # Get ligand name
    if re.search(r'^HETNAM', l):
      HetNam_Lns.append(l.rstrip())

####################
  ## Parse Mutation data found in PDB HEADER
  for l in SeqAdv:
    if re.search(r'MUTATION', l):
      Itms = l.split()
      chain, orig, posit, new = Itms[3], AA(Itms[7]), Itms[4], AA(Itms[2])
      change = '{0}:{1}{2}{3}'.format(chain, orig, posit, new)
      if mutation is None:
        mutation = ''.join(change)
      else:
        mutation += '|'+change

  ## Parse Heteratom data found in PDB HEADER: Ligand, unnatural amino acid
  Het_Uniq = list(set(Het_Lns))
  for het_id in Het_Uniq:
    het_nm = ''
    for l in HetNam_Lns:
      if re.search(het_id, l):
        het_nm += l.split(het_id)[1].strip()

    # format for ligand: <3-letter code>:<full name>
    lig_name = '{0}:{1}'.format(het_id, het_nm)
    # modified amino acid, save
    if len(ModRes) > 0:
      for mod_res in ModRes:
        if het_id == mod_res:
          if aa_modif is None:
            aa_modif = ''.join(lig_name)
          else:
            aa_modif += '|'+lig_name
          continue

    ## if it is salt ions and additives, etc
    if SaltAdditive(het_id):
      if UnnaturalAA(het_id):
        continue
      if salt is None:
        salt = ''.join(lig_name)
      else:
        salt += '|'+lig_name
    else:
      if ligand is None:
        ligand = ''.join(lig_name)
      else:
        ligand += '|'+lig_name

###############
  ## Identify individual protein chain from PDB HEADER, parse them individually
  ## and put their data into dictionary of object
  for l in headx:
    # Get and sort Protein/Organism data based on the MOL_ID, if multi-chain
    if re.search(r'^COMPND', l):
      # Set each found protein chain as individual object
      if re.search(r'MOL_ID:', l):
        mol_id = int(l.rstrip().split(':')[1].split(';')[0].strip())
        HD[mol_id] = HeadData(pmid=pmid, deposit=deposit, release=release,
                              latest=latest, resolu=resolu, space=space,
                              crystal=crystal, mutation=mutation,
                              ligand=ligand, salt=salt, aa_modif=aa_modif, 
                              seqadv=SeqAdv, dbref=DBRef)
      if re.search(r'MOLECULE:', l):
        HD[mol_id].p_name = l.rstrip().split(':')[1].split(';')[0].strip()
      if re.search(r'MUTATION:', l):
        HD[mol_id].mutate = l.rstrip().split(':')[1].split(';')[0].strip()
      if re.search(r'CHAIN:', l):
        HD[mol_id].chain = l.rstrip().split(':')[1].split(';')[0].strip()
#      if re.search(r'EC:', l):  EC found in PDB HEADER can be wrong in old PDBs
#        HD[mol_id].ec = l.rstrip().split(':')[1].split(';')[0].strip()

    # Get Gene data
    if re.search(r'^SOURCE', l):
      if re.search(r'MOL_ID:', l):
        mol_id = int(l.rstrip().split(':')[1].split(';')[0].strip())
      if re.search(r'ORGANISM_SCIENTIFIC:', l):
        HD[mol_id].species = l.rstrip().split(':')[1].split(';')[0].strip()
      if re.search(r'ORGANISM_COMMON:', l):
        HD[mol_id].common = l.rstrip().split(':')[1].split(';')[0].strip()
      if re.search(r'ORGANISM_TAXID:', l):
        HD[mol_id].taxid = l.rstrip().split(':')[1].strip(';').strip()
      if re.search(r'GENE:', l):
        HD[mol_id].gene = l.rstrip().split(':')[1].split(';')[0].strip()

####################################
  ## collect all info into datasets
  for chain_id in Entity.keys():
    for mol_id in HD.keys():
      if re.search(r'{0}'.format(Entity[chain_id].chain_id), HD[mol_id].chain):
        Entity[chain_id].p_name   = HD[mol_id].p_name
        Entity[chain_id].chain    = HD[mol_id].chain
        Entity[chain_id].mutate   = HD[mol_id].mutate
        Entity[chain_id].mutation = HD[mol_id].mutation
#        Entity[chain_id].ec       = HD[mol_id].ec
        Entity[chain_id].species  = HD[mol_id].species
        Entity[chain_id].common   = HD[mol_id].common
        Entity[chain_id].taxid    = HD[mol_id].taxid
        Entity[chain_id].gene     = HD[mol_id].gene
        Entity[chain_id].pmid     = HD[mol_id].pmid
        Entity[chain_id].deposit  = HD[mol_id].deposit
        Entity[chain_id].release  = HD[mol_id].release
        Entity[chain_id].latest   = HD[mol_id].latest
        Entity[chain_id].resolu   = HD[mol_id].resolu
        Entity[chain_id].space    = HD[mol_id].space
        Entity[chain_id].ligand   = HD[mol_id].ligand
        Entity[chain_id].salt     = HD[mol_id].salt
        Entity[chain_id].aa_modif = HD[mol_id].aa_modif
        Entity[chain_id].crystal  = HD[mol_id].crystal
        Entity[chain_id].uni_length = SearchUniProt( Entity[chain_id].uni_id )

  return Entity


##########################################################################
## Object to store Data. To make it a iterable dict in python3, output has iter()
class HeadData(object):
  def __init__(self,pdb_id=None, pdb=None, chain_id=None, conf=None, 
                    uni_id=None, pdb_length=None, uni_length=None,
                    p_name=None, chain=None, ec=None, gene=None,
                    species=None, common=None, taxid=None,
                    pmid=None, deposit=None, release=None, latest=None, 
                    resolu=None, mutate=None, mutation=None, space=None, 
                    ligand=None, salt=None, aa_modif=None,
                    crystal=None, seqadv=None, dbref=None ):
      self.pdb_id = pdb_id      # pdb_x ID for PDB and chain ID
      self.pdb = pdb            # pdb ID for PDB (orig input)
      self.chain_id = chain_id  # pdb chain ID for PDB (orig input)
      self.conf = conf          # xtal conformation from Konformation  
      self.uni_id = uni_id              # uniprot ID found in Entity query
      self.pdb_length = pdb_length      # chain length found in Entity query
      self.uni_length = uni_length      # (full) chain length found in uniprot
      self.p_name  = p_name     # protein name
      self.chain = chain        # chain ID in PDB
      self.mutate = mutate      # presence of mutation
      self.mutation = mutation  # mutation name
      self.ec    = ec           # EC number of protein
      self.species = species    # scientific name of organism of protein
      self.common = common      # common name of organism of protein
      self.taxid = taxid        # taximony ID
      self.gene  = gene         # Gene name
      self.pmid  = pmid         # PubMed ID
      self.resolu = resolu      # crystal resolution
      self.deposit = deposit    # original deposit date, not last update date
      self.release = release    # deposit release date, not last update date
      self.latest  = latest     # latest date, indicated update to deposit
      self.space = space        # xtal lattice space group
      self.ligand = ligand      # HETNAM in PDB header, organic chemicals
      self.salt = salt          # HETNAM in PDB header, salt ions 
      self.aa_modif = aa_modif  # MODRES in PDB header, unnatural amino acid
      self.crystal = crystal    # CRYSTL in PDB header, xtal lattice group info
      self.seqadv = seqadv      # SEQADV in PDB header, mutation info
      self.dbref = dbref        # DBREF in PDB header

  def __iter__(self):
    return iter( self.__dict__.items() )


##########################################################################
# Search UniProt webpage for full-length data; only use canonical info
def SearchUniProt( uni_id ):

  uni_length = None
  if uni_id is not None:
    html_uni = 'https://www.uniprot.org/uniprot/{0}.txt'.format(uni_id)
    uni_info   = requests.get(html_uni)

    ## downloaded FASTA is byte, need to split it into individual chains by '\n',
    unix = [ x for x in uni_info.content.decode().split('\n') if x is not '' ]
    for l in unix:
      if re.search('^SQ ', l):
        uni_length = int(l.split('SEQUENCE')[1].split('AA')[0].strip())

  return uni_length

###########################################################################
#
#  Peter M.U Ung @ MSSM/Yale
#
#  v1.0   20.01
#  v2.0   20.03.01  change result xml reading from string parsing to dict
#
#