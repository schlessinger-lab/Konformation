#!/usr/bin/python

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
from pathos import multiprocessing
from tqdm import tqdm
from aa_residue     import *
from x_data_coll    import *
from x_r_c_spines   import *
from x_dfg_torsion  import *
from x_helix_axis   import *
from x_domain_dist  import *
from x_pdb_extract  import *
from x_fasta_parse  import *
from x_ligand_type  import *
from x_konf_vars    import *
from x_search_align import *
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
p = PDBParser(PERMISSIVE=1, QUIET=True)

##########################################################################
def Konformation( RefRes, Query_PDB, output, **kwargs ):

  # Read in aligned fasta library for the PDB
  Fasta_Lib = CacheSeqDatabase(RefRes['FASTA'])
  print(' ** Found number of aligned FASTA entries: '+str(len(Fasta_Lib)))

  # Create master parameter database for all query PDBs, including Ref PDB
  PDB_Data = {}
  PDB_Data[RefRes['REFPDB'].split('/')[-1].split('.')[0]] = BuildDataSet()
  for pdb in Query_PDB:
    pdb_id = pdb.split('/')[-1].split('.')[0]
    PDB_Data[pdb_id] = BuildDataSet()
    PDB_Data[pdb_id]['pdb_id'] = pdb_id
    PDB_Data[pdb_id]['pdb']    = pdb


  ParameterCalculations( RefRes, Query_PDB, Fasta_Lib, PDB_Data, output )

  # Print out all data in dataframe
  Data2Pandas(PDB_Data, output)


##########################################################################
def ParameterCalculations( RefRes, Query_PDB, Fasta_Lib, PDB_Data, output ):

  ref_pdb   = RefRes['REFPDB']
  pdb_dir   = RefRes['PDBDIR']
  missing   = RefRes['MISSRES']
  helix_res = RefRes['HELIX']
  n_dom_res = RefRes['NDOM']
  c_dom_res = RefRes['CDOM']
  dfg_f_res = RefRes['DFGF']
  gate_res  = RefRes['GATE']

  xheli_res = RefRes['XHELIX']
  zndom_res = RefRes['ZNDOM']
  zcdom_res = RefRes['ZCDOM']
  zdfgf_res = RefRes['ZDFGF']

  # Read in Reference PDB, extract the helix residues and domain resid columns,
  ref_pdb_id = ref_pdb.split('/')[-1].split('.')[0]
  print('\n### Extract reference PDB coordinates and parameters ###')
  print('>>>> Reference PDB: '+ref_pdb_id)
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
    Ref_Gate[name]  = [ ref_gate[name], [] ]


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
  pRef = ParsePDB( h_seq=Ref_Helix, n_seq=Ref_N_Dom, c_seq=Ref_C_Dom,
                   g_seq=Ref_Gate,  f_seq=Ref_DFG_F, t_seq=Ref_Cs, 
                   pdb_dir=pdb_dir )
  Ref_Coords = pRef.extract_pdb(ref_pdb)
  
  # Calculate the query PDBs helix axis and C-helix parameters
  mpi  = multiprocessing.Pool(processes = multiprocessing.cpu_count())
  pPDB = ParsePDB( h_seq=Helix_Seq, n_seq=N_Dom_Seq, c_seq=C_Dom_Seq,
                   g_seq=Gate_Seq,  f_seq=DFG_F_Seq, t_seq=Cs_Seq,
                   pdb_dir=pdb_dir, corr=CoordCorrect(missing, pdb_dir) )
#  Tmp  = [pPDB(pdb) for pdb in Query_PDB]
  Tmp  = [x for x in tqdm(mpi.imap_unordered(pPDB,Query_PDB),total=len(Query_PDB))]
  mpi.close()
  mpi.join()

  # PDB_Coords = [pdb_name, H_Crds, N_Crds, C_Crds, G_Crds, F_Crds, T_Crds]
  PDB_Coords = [Itm for Itm in Tmp if Itm is not None]
#  PDB_Coords = Tmp
  print('\n## Input Query: {0} - Accepted: {1}\n'.format(
        len(Tmp), len(PDB_Coords)))
  os.system('cat _TEMP.missing.* > '+output+'.missing.txt; rm _TEMP.missing.*')

  # Compare reference and query PDB C-helix/N-dom/C-dom parameters
  RefReg2, Reg2 = HelixMeasurements(Ref_Coords, PDB_Coords, PDB_Data, output)

  DomainDistances( Ref_Coords, PDB_Coords, RefReg2, Reg2, PDB_Data, output)
  DFGTorsionAngle( Ref_Coords, PDB_Coords, PDB_Data, output)
  DescriptLigands( Ref_Coords, PDB_Coords, PDB_Data )


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
