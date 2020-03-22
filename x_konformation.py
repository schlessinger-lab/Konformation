#!/usr/bin/env python3

##########################################################################
##
##	Peter M.U. Ung	@ MSSM
##
##  ** Converted from the standalone '0_kinmetrics.py'
##
##	Calculate the helix axis of the kinase C-helix, and use the axis and
##	conserved Glu on C-helix to measure whether the kinase structure has
##	C-in or C-out conformation, relative to a reference inhibitor-bound 
##	C-in structure (1atp).
##	- C-helix axis angle
##	- C-helix curvature
##	- C-helix/N-domain Glu/Lys distance
##
## 	** Make sure no 'blank' residue in fasta library in the regions
##	   that will be used for residue extraction -- catalytic conserved 
##	   positions (+/- 3 residues) of C-helix (Glu) and N-domain (Lys)
##
##########################################################################

import sys,os
import re
from tqdm import tqdm
from pathos import multiprocessing

from x_pdb_extract  import ParsePDB
from x_data_coll    import Data2Pandas
from x_data_coll    import BuildDataSet
from x_pdb_extract  import CoordCorrect
from x_fasta_parse  import RefFastaColumn
from x_fasta_parse  import FastaFromColumn
from x_ligand_type  import DescriptLigands
from x_dfg_torsion  import DFGTorsionAngle
from x_domain_dist  import DomainDistances
from x_search_align import CacheSeqDatabase
from x_helix_axis   import HelixMeasurements
#from x_x_r_c_spines import RCSpinesMeasure

from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
p = PDBParser(PERMISSIVE=1, QUIET=True)

##########################################################################
def Konformation( parm, Query_PDB, output, **kwargs ):

  # Read in aligned fasta library for the PDB
  Fasta_Lib = CacheSeqDatabase(parm['FASTA'][0])
  print(' ** Found number of aligned FASTA entries: \033[31m{0}\033[0m'.format(len(Fasta_Lib)))

  # Create master parameter database for all query PDBs, including Ref PDB
  PDB_Data = {}
  PDB_Data[parm['REFPDB'][0].split('/')[-1].split('.')[0]] = BuildDataSet()
  for pdb in Query_PDB:
    pdb_id = pdb.split('/')[-1].split('.')[0]
    PDB_Data[pdb_id] = BuildDataSet()
    PDB_Data[pdb_id]['pdb']    = pdb
    PDB_Data[pdb_id]['pdb_id'] = pdb_id

  ## Calculated data is stored back into the same PDB_Data
  ParameterCalculations( parm, Query_PDB, Fasta_Lib, PDB_Data, output )

  # Print out all data in dataframe
  parm_df = Data2Pandas(PDB_Data, output)
  
  return parm_df


##########################################################################
## Calculated data is stored back into the same PDB_Data
def ParameterCalculations( parm, Query_PDB, Fasta_Lib, PDB_Data, output ):

  ref_pdb   = parm['REFPDB'][0]
  pdb_dir   = parm['PDBDIR'][0]
  missing   = parm['MISSRES'][0]
  ## primary seq sets for key recognition sites
  helix_res = parm['HELIX']
  n_dom_res = parm['NDOM']
  c_dom_res = parm['CDOM']
  dfg_f_res = parm['DFGF']
  gate_res  = parm['GATE']

  ## if failed with primary set, try alternative seq that takes last residues
  ## alternative seq sets for key recognition sites, the last resid
  xheli_res = parm['XHELIX']
  zndom_res = parm['ZNDOM']
  zcdom_res = parm['ZCDOM']
  zdfgf_res = parm['ZDFGF']

  # Read in Reference PDB, extract the helix residues and domain resid columns,
  ref_pdb_id = ref_pdb.split('/')[-1].split('.')[0]
  print('\n\033[34m### Extract reference PDB coordinates and parameters ###\033[0m')
  print('>>>> Reference PDB: \033[31m{}\033[0m'.format(ref_pdb_id))
  ref_helix, helix_column = RefFastaColumn(Fasta_Lib, helix_res, ref_pdb_id)
  ref_n_dom, n_dom_column = RefFastaColumn(Fasta_Lib, n_dom_res, ref_pdb_id)
  ref_c_dom, c_dom_column = RefFastaColumn(Fasta_Lib, c_dom_res, ref_pdb_id)
  ref_dfg_f, dfg_f_column = RefFastaColumn(Fasta_Lib, dfg_f_res, ref_pdb_id)
  ref_gate,  gate_column  = RefFastaColumn(Fasta_Lib, gate_res,  ref_pdb_id)

  xref_heli, xheli_column = RefFastaColumn(Fasta_Lib, xheli_res, ref_pdb_id)
  zref_ndom, zndom_column = RefFastaColumn(Fasta_Lib, zndom_res, ref_pdb_id)
  zref_cdom, zcdom_column = RefFastaColumn(Fasta_Lib, zcdom_res, ref_pdb_id)
  zref_dfgf, zdfgf_column = RefFastaColumn(Fasta_Lib, zdfgf_res, ref_pdb_id)

  Ref_Helix, Ref_N_Dom, Ref_C_Dom, Ref_DFG_F, Ref_Gate = {}, {}, {}, {}, {}
  for name in ref_helix:
    Ref_Helix[name] = [ ref_helix[name], xref_heli[name] ]
    Ref_N_Dom[name] = [ ref_n_dom[name], zref_ndom[name] ]
    Ref_C_Dom[name] = [ ref_c_dom[name], zref_cdom[name] ]
    Ref_DFG_F[name] = [ ref_dfg_f[name], zref_dfgf[name] ]
    Ref_Gate[name]  = [ ref_gate[name],  [] ]


  # Extract residues using the column info
  helix_seq = FastaFromColumn(Fasta_Lib, helix_column, len(helix_res))
  n_dom_seq = FastaFromColumn(Fasta_Lib, n_dom_column, len(n_dom_res))
  c_dom_seq = FastaFromColumn(Fasta_Lib, c_dom_column, len(c_dom_res))
  dfg_f_seq = FastaFromColumn(Fasta_Lib, dfg_f_column, len(dfg_f_res))
  gate_seq  = FastaFromColumn(Fasta_Lib, gate_column,  len(gate_res) )

  xHeli_seq = FastaFromColumn(Fasta_Lib, xheli_column, len(xheli_res))
  zNDom_seq = FastaFromColumn(Fasta_Lib, zndom_column, len(zndom_res))
  zCDom_seq = FastaFromColumn(Fasta_Lib, zcdom_column, len(zcdom_res))
  zDFGF_seq = FastaFromColumn(Fasta_Lib, zdfgf_column, len(zdfgf_res))

  Helix_Seq, N_Dom_Seq, C_Dom_Seq, DFG_F_Seq, Gate_Seq = {}, {}, {}, {}, {}
  for name in helix_seq:
    Helix_Seq[name] = [ helix_seq[name], xHeli_seq[name] ]
    N_Dom_Seq[name] = [ n_dom_seq[name], zNDom_seq[name] ]
    C_Dom_Seq[name] = [ c_dom_seq[name], zCDom_seq[name] ]
    DFG_F_Seq[name] = [ dfg_f_seq[name], zDFGF_seq[name] ]
    Gate_Seq[name]  = [ gate_seq[name], [] ]

  # Extract R/C-spine residues
#  Ref_Rs, Rs_Seq  = [], []
  Ref_Cs, Cs_Seq  = [], []


  # Calculate the reference helix axis and C-helix parameters
  pRef = ParsePDB(  h_seq=Ref_Helix, n_seq=Ref_N_Dom, c_seq=Ref_C_Dom,
                    g_seq=Ref_Gate,  f_seq=Ref_DFG_F, t_seq=Ref_Cs, 
                    pdb_dir=pdb_dir, corr='' )
  Ref_Coords = pRef.extract_pdb(ref_pdb)
  
  # Calculate the query PDBs helix axis and C-helix parameters
  pPDB = ParsePDB(  h_seq=Helix_Seq, n_seq=N_Dom_Seq, c_seq=C_Dom_Seq,
                    g_seq=Gate_Seq,  f_seq=DFG_F_Seq, t_seq=Cs_Seq,
                    pdb_dir=pdb_dir, corr=CoordCorrect(missing, pdb_dir) )

  # Extract PDB data. Using MPI can slow down alot due to IO limit
  Tmp  = [ pPDB(pdb) for pdb in tqdm(Query_PDB) ]

  # PDB_Coords = [pdb_name, H_Crds, N_Crds, C_Crds, G_Crds, F_Crds, T_Crds]
  PDB_Coords = [Itm for Itm in Tmp if Itm is not None]

  print('\n ## Input Query: {0} - Accepted: \033[31m{1}\033[0m\n'.format(
        len(Tmp), len(PDB_Coords)))

  ## Handle any structure that is flagged with missing data
  missing = False
  for File in os.listdir('.'):
    if re.search('_TEMP.missing.', File):
      missing = True
      with open('_TEMP.missing-header.txt', 'w') as fo:
        fo.write('## Recover missing residue coordinates if *middle* residue is available\n')
        fo.write('## by saving residue coordinates in PDB format: correct.<domain>.<full_pdb_name>.pdb\n')
        fo.write('## then rerun the script with modified input with updated MISSRES <new list filename>\n')
        fo.write('##PDB_File|Missing_domain|residues\n')
  if missing:
    print('\033[31m## Check for missing structural data:\033[0m '+output+'.missing.txt')
    os.system('cat _TEMP.missing-header.txt _TEMP.missing.* > '+output+'.missing.txt; rm _TEMP.missing.* _TEMP.missing-header*')


  # Compare reference and query PDB C-helix/N-dom/C-dom parameters
  RefReg2, Reg2 = HelixMeasurements( Ref_Coords, PDB_Coords, PDB_Data, parm, output )

  ## Measeure various parameters; RC-Spines are irrelevant and are not calculated anymore
  DomainDistances( Ref_Coords, PDB_Coords, RefReg2, Reg2, PDB_Data, parm, output )
  DFGTorsionAngle( Ref_Coords, PDB_Coords, PDB_Data, parm, output )
  DescriptLigands( Ref_Coords, PDB_Coords, PDB_Data )
#  RCSpinesMeasure( Ref_Coords, PDB_Coords, PDB_Data, parm, output )

##########################################################################
##
##  v0.0    17.01.08
##  v0.1    17.01.23    moved sub functions to separate files
##  v0.2    17.01.28
##  v0.3    17.01.28
##
##  v1.0    17.02.02
##  v1.1    17.02.02	minor changes
##  v2.0    17.02.22	add data storage object to collect data
##  v3.0    17.03.08    modified to include R/C-spine calculation
##  v4.0    17.03.24    simplified input with use of parameter file
##  v5.0    17.04.12    moved optional missing file to parameter file
##  v6.0    17.08.18	add ligand type measurement
##  v7.0    17.12.21    only calculate essential parameters
##  v8.0    18.03.12    change it into a callable and accept PDB w/ full path
##  v9.0    18.03.15    masked all R/C-spine calculation, improve in all calc
##
##########################################################################
