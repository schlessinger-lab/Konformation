#!/usr/bin/env python3

import re
import os
import time
import requests

import numpy as np
import pandas as pd

from aa_residue import AA
from aa_residue import UnnaturalAA
from aa_residue import SaltAdditive

from tqdm import tqdm

##########################################################################   
## Read in pdb list with different formats:
# Format of the list:
#  1) <pdb_id>_<chain_id>.xxx.pdb
#  2) <pdb_id>_<chain_id>.xxx.pdb <chain_id>
#  3) <pdb_id>.xxx.pdb            <chain_id>
#  4) <pdb_id>                    <chain_id>
#  (PDB ID and chain ID separated by '_')
#
#  e.g.:   1ATP_E.xxx.pdb
#          3HHP.xxx.pdb    C
#          6GTT            A
def ProcessPDBList( pdb_list ):
  PDB_List = pd.read_csv(pdb_list, sep='\s+', header=None, comment='#').to_numpy()

  PDB   = []
  for entry in PDB_List:
    if len(entry) == 2:
      if re.search(r'.pdb', entry[0]):
        pdb_id   = entry[0].split('.')[0].split('_')[0]
        chain_id = entry[0].split('.')[0].split('_')[1]
      else:
        pdb_id   = entry[0]
      if type(entry[1]) == str:
        chain_id = entry[1]
      PDB.append([pdb_id, chain_id])
    elif len(entry) == 1:
      if re.search(r'.pdb', entry):
        pdb_id   = entry.split('.')[0].split('_')[0]
        chain_id = entry.split('.')[0].split('_')[1]
      else:
        pdb_id   = entry.split('_')[0]
        chain_id = entry.split('_')[1]
      PDB.append([pdb_id, chain_id])
    else:
      sys.exit('\033[31m ERROR: PDB list format is invalid\033[0m')

  return PDB


##########################################################################
## Download the PDB page from the Web and search for the UniProt ID
def SearchPDB( pdb ):
  pdb_id, chain_id = pdb
#  print(pdb_id, chain_id)

  pdb_length, uni_id, tax_id, species = None, None, None, None
  ec, pmid, p_name, uni_id = None, None, None, None
  resolu, deposit, release, latest = None, None, None, None

  ## Retreive basic PDB data from RCSB query search
  ## delay retrieve by 0.1s to avoid server from download lockout (prevent downloading)
  html_pdb   = 'https://www.rcsb.org/pdb/rest/describePDB?structureId={0}'.format(pdb_id)
  html_chain = 'https://www.rcsb.org/pdb/rest/describeMol?structureId={0}.{1}'.format(pdb_id, chain_id)
  pdb_info   = requests.get(html_pdb)
  time.sleep(0.1)
  chain_info = requests.get(html_chain)

  ## downloaded FASTA is byte, need to split it into individual chains by '>',
  ## then concatenate pdb_id and chain_id. Split string by 'SEQUENCE', then
  ## remove all '\n'
  pdbx = [ re.sub('/|<','', x).strip() for x in pdb_info.content.decode().split('>') if x is not '' ]
  chnx = [ re.sub('/|<','', x).strip() for x in chain_info.content.decode().split('>') if x is not '' ]

  for l in pdbx:
    if re.search('structureId', l):
      if re.search('pubmedId=', l):
        pmid   = re.sub('"', '', l.split('pubmedId=')[1].split()[0])
      if re.search('resolution=', l):
        resolu = re.sub('"', '', l.split('resolution=')[1].split()[0])
      if re.search('deposition_date=', l):
        deposit= re.sub('"', '', l.split('deposition_date=')[1].split()[0])
      if re.search('release_date=', l):
        release= re.sub('"', '', l.split('release_date=')[1].split()[0])
      if re.search('last_modification_date=', l):
        latest = re.sub('"', '', l.split('last_modification_date=')[1].split()[0])

  for l in chnx:
    if re.search('length=', l):
      pdb_length = int(re.sub('"', '', l.split('length=')[1].split()[0]))
    if re.search('Taxonomy ', l):
      tax_id  = re.sub('"', '', l.split(' id=')[1])
      species = re.sub('"', '', l.split(' name=')[1].split(' id')[0])
    if re.search('macroMolecule ', l):
      p_name  = re.sub('"', '', l.split(' name=')[1])
    if re.search('accession id=', l):
      uni_id  = re.sub('"', '', l.split(' id=')[1])
    if re.search('enzClass ec=', l):
      ec = re.sub('"', '', l.split(' ec=')[1])

#########################################
  ## Find UniProt ID, chain ID, chain Length, in PDB webpage

  Entity = {}
  Entity[chain_id] = HeadData(pdb_id=pdb_id, chain_id=chain_id, ec=ec,
                              pdb_length=pdb_length, uni_id=uni_id,
                              taxid=tax_id, species=species, p_name=p_name,
                              pmid=pmid, resolu=resolu, latest=latest,
                              deposit=deposit, release=release)
  Entity[chain_id].chain = chain_id

################################
  # Retrieve data from HEADER data in RCSB webpages
  html_head = 'https://files.rcsb.org/header/{0}.pdb'.format(pdb_id)
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

    ## if it is salt ions and additives, etc
    if SaltAdditive(het_id):
      if salt is None:
        salt = ''.join(lig_name)
      else:
        salt += '|'+lig_name
    ## if it is unnatural Amino Acids, ignore -- registered by aa_modif
    elif UnnaturalAA(het_id):
      continue    # if ligands, save; multiple ligands, use '|' as separator
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
        Entity[chain_id].uni_length = SearchUniProt( Entity[chain_id].uni_id,
                                                Entity[chain_id].pdb_id )

  return Entity


##########################################################################
## Object to store Data. To make it a iterable dict in python3, output has iter()
class HeadData(object):
  def __init__(self, pdb_id=None, chain_id=None, uni_id=None,
                    pdb_length=None, uni_length=None,
                    p_name=None, chain=None, ec=None, gene=None,
                    species=None, common=None, taxid=None,
                    pmid=None, deposit=None, release=None, latest=None, 
                    resolu=None, mutate=None, mutation=None, space=None, 
                    ligand=None, salt=None, aa_modif=None,
                    crystal=None, seqadv=None, dbref=None ):
      self.pdb_id = pdb_id      # pdb ID for PDB (orig input)
      self.chain_id = chain_id    # pdb chain ID for PDB (orig input)
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
def SearchUniProt(uni_id, pdb_id):

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