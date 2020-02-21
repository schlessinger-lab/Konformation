#!/usr/bin/env python3

import sys

##########################################################################
##
##	Peter M.U. Ung	@ MSSM
##
##  Deceptrated: superseded by 1_run_konformation.py
##
##  **  '0_kinmetrics-orig.py' - the original script to run descriptor
##                               generation with all kinase PDB in a folder
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
from x_konf_vars import GenerateTemplSetupScript

templ_file = 'Template_parameter_file.in'
msg = '''\n    Usage: {0}             
            [Parameter file]
	    [Output filename]

            ( correction file <-- corrected file: correct.<PDB file> )
                                  Will need to rerun 2-3x to confirm 
                                  no NEW missing residue
            (2-3 residues before and after the Center residue (Glu/Lys/Asp))\n
   e.g.: > {0}\n              parameter.file output
\t ** Generated - \033[31m{1}\033[0m\n'''.format(sys.argv[0], templ_file)
if len(sys.argv) != 3 :
  GenerateTemplSetupScript(templ_file)
  sys.exit(msg) 

##########################################################################
import re,os,glob
from tqdm import tqdm
from pathos import multiprocessing

from x_data_coll   import BuildDataSet
from x_data_coll   import Data2Pandas
#from x_r_c_spines  import *
from x_dfg_torsion import DFGTorsionAngle
from x_helix_axis  import HelixMeasurements
from x_domain_dist import DomainDistances
from x_pdb_extract import ParsePDB
from x_pdb_extract import CoordCorrect
from x_fasta_parse import RefFastaColumn
from x_fasta_parse import FastaFromColumn
from x_ligand_type import DescriptLigands
from x_pdb_superpose import SuperposePDB
from x_konformation  import Konformation
from x_konf_vars     import ParseParameterFile
from x_search_align  import CacheSeqDatabase
from x_search_align  import CheckInputStructures
from x_search_align  import GenerateProfileAlignment

from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
p = PDBParser(PERMISSIVE=1, QUIET=True)

##########################################################################
def main( param_list, output, **kwargs ):

  # Read in parameter files for residues and other params
  RefRes = ParseParameterFile(param_list)

  # Read in aligned fasta library for the PDB as dictionary
  Fasta_Lib = {}
  print('\n### Reading in aligned FASTA library ###')
  Temp_Lib = SeqIO.to_dict(SeqIO.parse(RefRes['FASTA'], 'fasta'))
  print(' Found aligned FASTA entries:\t'+str(len(Temp_Lib)))
  for key in Temp_Lib.keys():
    new_key = key.split('|')[0]
    Fasta_Lib[new_key] = Temp_Lib[key]


  ref_pdb = RefRes['REFPDB']

  # Read in Query PDB names from the input list
  with open(RefRes['PDBLIST'], 'r') as fi:
    Query_PDB = [line.rsplit()[0] for line in fi]

  # Create the master database for parameters of all query PDBs
  PDB_Data = {}
  PDB_Data[ref_pdb.split('.')[0]] = BuildDataSet()
  for pdb in Query_PDB:
    pdb_id = pdb.split('.')[0]
    PDB_Data[pdb_id] = BuildDataSet()
    PDB_Data[pdb_id]['pdb_id'] = pdb_id


  ParameterCalculations( Fasta_Lib, Query_PDB, RefRes, PDB_Data,
                         ref_pdb, output )

  # Print out all data in dataframe
  Data2Pandas(PDB_Data, output)


##########################################################################
##########################################################################
def ParameterCalculations( Fasta_Lib, Query_PDB, RefRes, PDB_Data, 
          ref_pdb, output, **kwargs):

  pdb_dir   = RefRes['PDBDIR']
  missing   = RefRes['MISSRES']
  helix_res = RefRes['HELIX']
  n_dom_res = RefRes['NDOM']
  c_dom_res = RefRes['CDOM']
  gate_res  = RefRes['GATE']
  Rs_List   = []
  Cs_List   = []
#  Rs_List = [ RefRes['RSPINE1'], RefRes['RSPINE2'],
#              RefRes['RSPINE3'], RefRes['RSPINE4'] ]
#  Cs_List = [ RefRes['CSPINE1'], RefRes['CSPINE2'], RefRes['CSPINE3'],
#              RefRes['CSPINE4'], RefRes['CSPINE5'], RefRes['CSPINE6'],
#              RefRes['CSPINE7'], ]# RefRes['CSPINE8'] ]

  # Read in Reference PDB, extract the helix residues and domain resid columns,
  ref_pdb_name = ref_pdb.split('.')[0]
  print('\n### Extract reference PDB coordinates and parameters ###')
  print('>>>> Reference PDB: '+ref_pdb_name)

  ## Extract sequence in reference PDB using column info: dict of seq
  Ref_Helix, helix_column = RefFastaColumn(Fasta_Lib, helix_res, ref_pdb_name)
  Ref_N_Dom, n_dom_column = RefFastaColumn(Fasta_Lib, n_dom_res, ref_pdb_name)
  Ref_C_Dom, c_dom_column = RefFastaColumn(Fasta_Lib, c_dom_res, ref_pdb_name)
  Ref_Gate,  gate_column  = RefFastaColumn(Fasta_Lib, gate_res,  ref_pdb_name)

  ## Extract sequence in unknown PDB using column info: dict of seq
  Helix_Seq = FastaFromColumn(Fasta_Lib, helix_column, len(helix_res))
  N_Dom_Seq = FastaFromColumn(Fasta_Lib, n_dom_column, len(n_dom_res))
  C_Dom_Seq = FastaFromColumn(Fasta_Lib, c_dom_column, len(c_dom_res))
  Gate_Seq  = FastaFromColumn(Fasta_Lib, gate_column,  len(gate_res) )

  ## R- and C-spines, but useless so skip them: list(dict of seq)
  Ref_Rs, Rs_Seq  = [], []
  Ref_Cs, Cs_Seq  = [], []
  for idx, residues in enumerate(Rs_List):   # [Ref_Residues, Ref_column]
    Rs_Tmp, col = RefFastaColumn(Fasta_Lib, residues, ref_pdb_name)
    Ref_Rs.append(Rs_Tmp)
    Rs_Seq.append(FastaFromColumn(Fasta_Lib, col, len(Rs_List[idx])))

  for idx, residues in enumerate(Cs_List):
    Cs_Tmp, col = RefFastaColumn(Fasta_Lib, residues, ref_pdb_name)
    Ref_Cs.append(Cs_Tmp)
    Cs_Seq.append(FastaFromColumn(Fasta_Lib, col, len(Cs_List[idx])))

  ## Calculate the reference helix axis and C-helix parameters
  pRef = ParsePDB( h_seq=Ref_Helix, n_seq=Ref_N_Dom, c_seq=Ref_C_Dom,
                   g_seq=Ref_Gate,  f_seq=Ref_Rs,    t_seq=Ref_Cs,
                   pdb_dir=pdb_dir )
  Ref_Coords = pRef.extract_pdb(ref_pdb)
  
  ## Calculate the query PDBs helix axis and C-helix parameters
  mpi  = multiprocessing.Pool()
  pPDB = ParsePDB( h_seq=Helix_Seq, n_seq=N_Dom_Seq, c_seq=C_Dom_Seq,
                   g_seq=Gate_Seq,  f_seq=Rs_Seq,    t_seq=Cs_Seq,
                   pdb_dir=pdb_dir, corr=CoordCorrect(missing, pdb_dir) )
  Tmp  = [pPDB(pdb) for pdb in Query_PDB]
#  Tmp  = [x for x in tqdm(mpi.imap(pPDB, Query_PDB), total=len(Query_PDB))]
  mpi.close()
  mpi.join()

  # PDB_Coords = [pdb_name, H_Crds, N_Crds, C_Crds, G_Crds, R_Crds, T_Crds]
  PDB_Coords = [Itm for Itm in Tmp if Itm is not None]
  print(PDB_Coords)
#  PDB_Coords = Tmp
  print('\n## Input Query: {0} - Accepted: {1}\n'.format(
        len(Tmp), len(PDB_Coords)))
  os.system('cat _TEMP.* > '+output+'.missing.txt; rm _TEMP.*')

  # Compare reference and query PDB C-helix/N-dom/C-dom parameters
  RefReg2, Reg2 = HelixMeasurements(Ref_Coords, PDB_Coords, PDB_Data, output)

  DomainDistances( Ref_Coords, PDB_Coords, RefReg2, Reg2, PDB_Data, output)
  DFGTorsionAngle( Ref_Coords, PDB_Coords, PDB_Data, output)
#  RCSpinesMeasure( Ref_Coords, PDB_Coords, PDB_Data, output)  
  DescriptLigands( Ref_Coords, PDB_Coords, PDB_Data ) 


##########################################################################
if __name__ == '__main__':
  main( sys.argv[1], sys.argv[2])

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
##
##########################################################################
