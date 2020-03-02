# Kinformation
**Classification of Protein Kinase (Catalytic Domain) Conformations**


__Citation 1:__ [\*Ung PMU, \*Rahman R, Schlessinger A. Redefining the Protein Kinase Conformational Space with Machine Learning. Cell Chemical Biology (2018) 25(7), 916-924.](https://doi.org/10.1016/j.chembiol.2018.05.002)

__Citation 2:__ [\*Rahman R, \*Ung PMU, Schlessinger A. KinaMetrix: a web resource to investigate kinase conformations and inhibitor space. 
Nucleic Acids Research (2019) 47(D1), D361–D366.](https://doi.org/10.1093/nar/gky916)

Reference for Modi-Dunbrack alignment: [Modi V, Dunbrack RL. A structurally-validated multiple sequence alignment of 497 human protein kinase domains. Scientific Reports (2019) 9, 19790.](https://doi.org/10.1038/s41598-019-56499-4)

Reference for Muscle MSA alignment: [MUSCLE: a multiple sequence alignment method with reduced time and space complexity. BMC Bioinf. (2004) 5, 113.](https://doi.org/10.1186/1471-2105-5-113)

Reference for BioPython: [Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics (2009) 25, 1422-3.](https://doi.org/10.1093/bioinformatics/btp163)

#########################################################################
- **Update the PDB structure library of Newly Published Protein Kinase Catalytic Domain and Retrieve Them**
```
1_find_new_kinase_pdb.py
    -r <setup name>  [ Read setup file, search RCSB PDB for new kinases ]

Optional:
    -set             [ Generate a default setup file used by -r ]
    -p <list name>   [ Skip RCSB search, download PDB of kinase listed ]

 e.g.> *.py -r setup.file
```
- This script searches the RCSB Protein Data Bank for all kinase-associated structures, then compare the query results (PDB_ID) to an existing internal kinase PDB library (known_list) to find those that aren't in the library yet. Since RCSB search query is 'dirty' and includes anything that has association to kinase (kinase subdmain, or has the name "kinase" in the Abstract), these new entries are checked by their FASTA sequence (default seq-length > 220 + seq-identity > 40% to human kinome) before downloading the kinase structures. The downloaded structures are then superposed onto 1ATP using a standard reference residue set.

- The new results are added to the existing internal kinase library list (known_list) for future reference. The script will not check and download structures with PDB and chain ID that have been checked before to avoid redundent download.

- As long as the (known_list) is kept current and the downloaded PDBs are placed into the correct folder, this script can be made into a repeating "cron" job to update the kinase structure library. **However**, the resultant **sequence alignment** and **superposed structures** should be **manually checked** to confirm there is no error.

- I tried using MPI to speed up the download of FASTA and PDB but seems to run into server lockdown if multiple jobs are retrieving from the same IP. Could be a defense mechanism against DOS attack?? Using serial retrieve with a delay of 0.1s in between seem to work fine. It took ~ 10 mins to finish the checking and downloading of hundreds of new structures after 2.5 yrs of no update. But if the update is kept up every few months, the process should be much faster.

- Multi-sequence alignment settings have been modified so that most new sequences aligned to a pre-aligned sequences (with gaps in precise positions) should be confined to a specifc spacing and gapping, where key positions, i.e. C-helix, DFG- and APE-motifs, are all aligned to the exact columns as the reference (1ATP) and kinome sequences. This is mostly acheived by changing the penalty points to gap opening and extension during the MSA profile alignment ("muscle -gapopen -5.0 -gapextend -2.0 -center 0.0"). So far all alignments have the same spacing and positions. With this setting, it is now possible to align new structures' sequence to the canonical human kinome alignment (originally used my own modified alignment, now switching to Modi-Dunbrack's), which results in alignment that enables automated extraction of key residue's adjacent sequences for use in conformation determination (see below). However, if the new kinase structures have insertion, the MSA might still be thrown off.

##########################################################################
- **Calculate kinase structural parameters from structures with pre-aligned sequences, Followed by ML-based Conformation Classification**
```
2_kinase_conf_classifier.py
    <file>   [ Parameter File ]

Optional:
    -t       [ Generate template input file ]

 e.g.: *.py parameter.file
       *.py -t     -> generates "Template_parameter_file.in"

Format of PDB list, as as "WTNEWKINASE" setting from 1_find_new_kinase_pdb.py:
  1) <pdb_id>_<chain_id>.xxx.pdb
```

- This is the main script to generate kinase structural parameters that would be used to generate classification of kinase structure. From a list of kinase structures (e.g. retreived from **1_find_new_kinase_pdb.py**, written out as "WTNEWKINASE" setting), there are 2 options to handle the PDB -- 1) as pre-superposed structures with pre-aligned sequences, 2) as un-superposed structures with un- or pre-aligned sequences. This script will align and superose the input PDB as needed, although this is _not recommended_ since there is a slight chance that the multi-seq alignment or the superposition might be off and uncorrected.

- Sometimes kinase PDB structures have missing electron density and thus missing residues that may be used by this script to determine the position of the key residues. This script will skip and flag down these problematic structures in the **first pass** in generate a file (with extension **<OUTPREF>.missing.txt**) with bad PDB and what issue it has. User will need to use the next script, **3_make_correction_pymol.py** to check and generate the coordinates of the missing residues. Once collected these missing coordinates, the user can run a **second pass** with an updated parameter file to indicate the list of these missing residue files (these missing residues PDBs should be next to the downloaded PDB).

- This script extracts the coordinates of key residues and calculate the structural parameters, which are then used in the conformation classification. This script generates **two** result files, one with structural parameters and another one with conformation of the input PDB structures.

- The _default_ Machine Learning algorithm used here is the **SKLearn** version of the **Random Forest (RF)** model, instead of the _R RF_ model used in the original article. However, since both R and SKLearn RF models are trained on the same dataset and used similar settings, the results from both models are very similar and should be interchangable. 
- In addition to _SKLearn RF_ model, I have created additional models that have similar or better accuracy to the original R RF model: _Neural Network (nn)_, _Gradient Boosting (gb)_, _K-nearest Neighbors (kn)_, _Decision Tree (dt)_, _Gaussian Process (gp)_, and _Support-Vector Machine (svm)_ (see my [/Kinformation_MD GitHub](https://github.com/mungpeter/Kinformation_MD) for out-of-bag errors). It appears _nn, kn, gb, svm_ models perform a slight bit better than _rf, dt_ models.

- Note, the current SKL ML models used here are from SKL _0.22.1_, which may not be backward compatible with newer SKL version in the future. In such even, will need to use [/Kinformation_MD/0_kinfo_SK_model_gen.py](https://github.com/mungpeter/Kinformation_MD/0_kinfo_SK_model_gen.py) to generate newer version of MK model for use.

- R randomforest algorithm is now removed since the package **rpy2** is not well maintained.

########
- Sometimes when there are many entries to run (like ~ 1000), some of the good structures might fail for no obvious reason (e.g. 6G36, 6G39, 6GTT), but they are fine when run with just a few entries (< 100).
```
4_Konformation/x_helix_axis.py:174: RankWarning: Polyfit may be poorly conditioned
  Fn2 = [LsqFit(list(range(xcount)), Coords[m:m-posit, x],2) for x in range(3)]

  #2# Domain Warning: Reg2 (C-helix linear regression) failed: 6G36_A
```
- This error happens in the Numpy **Polyfit** linear regression function, and my thought is it is a bug that has something to do with memory?? Since this issue pops up randomly, you will just have to watch for this error and redo the structures that got caught by this bug...


##########################################################################
- **Manually extact the missing residues from kinase PDB**
```
3_make_correction_pymol.py
  [ List of PDB with missing structural info ]

 e.g.>  x.py   <OUTPREF>.missing.txt
```

- This script creates a PyMOL session file with kinase structures that **2_kinase_conf_classifier.py** could not get the required structural data (either because of misalignment, failed superposition, missing residues, unrecognizable residue pattern).
- From the PyMOL, user will inspect the structures, manually select the residues required for structural parameter generation and then **rerun** *2_kinase_conf_classifier.py* to incorporate the missing residues.
- To get the residues, in the PyMOL session, select the residues as indicated by the type of missing info (GATE (gatekeeper), N_Dom (b3-Lys), C_Dom (DFG-D), DFG_F (DFG-F), CHELIX) and the flanking residues (2 residues on each side; C-helix Glu is better to have up to 4 residues on each side). Save the residues by the PyMOL command:
```
  pymol> save correct.<missing_info>.<pdb_file>, sele
e.g.)
  pymol> save correct.GATE.5GGT_A.1atp.pdb, sele
```
- All resultant correction PDB files should be made into a list and supply to the input parameter file's **MISSRES** flag, and supply **MISSREF** flag with the file _<OUTPREF>.missing.txt_ to indicate what have been corrected.

##########################################################################
- **Classification of kinase conformation with ML Algorithm from ONLY the pre-calculated structural parameters**
```
4_kinase_conf_ML_only.py
    -data   <file>      [ CSV file of re-calculated PDB structural parameters ]
    -mdldir <path>      [ Path to directory with ML models ]
    -out    <prefix>    [ Output prefix ]

  Optional:
    -use_sk <model>     [ Use SKLearn ML model: rf|svm|nn|kn|dt|gp|gb (def: rf) ]

e.g.>  4_kinase_conf_ML_only.py
          -data data.csv  -use_sk gb
          -mdldir '/Users/xxx/scripts/Kinformation/z_database'
```

- This is a strip-down version of **2_kinase_conf_classifier.py** as it **only** does classification of kinase conformation from a csv with pre-calculated kinase structural parameters.
- R randomforest algorithm is now removed since the package **rpy2** is not well maintained.


##########################################################################
- **Extract PDB Header and Ligand information**
```
5_extract_pdb_header_info.py
    [ list of PDB ]
    [ output prefix ]

Format of the list:
  1) <pdb_id>_<chain_id>.xxx.pdb
  2) <pdb_id>_<chain_id>.xxx.pdb <chain_id>
  3) <pdb_id>.xxx.pdb            <chain_id>
  4) <pdb_id>                    <chain_id>
  (PDB ID and chain ID separated by '_')

  e.g.:   1ATP_E.xxx.pdb
          3HHP.xxx.pdb    C
          6GTT            A

 e.g.> 5_extract_pdb_header_info.py 
        pdb.list   
        pdb_info

columns = [ 'pdb_id',  'chain_id', 'pdb_length', 'uni_length', 'uni_id', 
            'gene',    'p_name',   'mutate',     'mutation',   'ec',
            'species', 'common',   'taxid',      'deposit',    'release', 'latest',
            'ligand',  'salt',     'aa_modif',   'resolu',     'space',   'pmid' ]
```

- From a list of PDB files (with chain ID) to extract structure information, ligand name, etc. Search the PDB chain length and the corresponding full-length from RCSB PDB and UniProt databases via internet. No need to read from hard-copy of PDB files.

- The search process cannot use MPI, RCSB blocks multiple requests from same IP address simultaneous? (prevent DoS attack?) 

##########

```
6_extract_pdb_ligand.py
    [ Input file of Parsed PBB Header file with "pdb_id" and "ligand" ]
    [ Output prefix ]

 e.g.> ./6_extract_pdb_ligand.py
          test-sample.example.csv.gz  test-again
```

- From an input file of Parsed PDB Header file with information on  "pdb_id" and "ligand" (result from **5_extract_pdb_header_info.py**), download the SMILES string of bound ligands and associate these ligands to the "pdb_id" they came from.

##########################################################################
- **Example 1: Retreive new kinase structures that is not in PDB collection**
```
> examples/
    |--1_find_new_kinase_pdb
            |------ run_find_new_kinase.sh    # running this setup
            |------ check_kinase.200227.fasta               # all aligned kinase PDB seq
            |------ check_kinase.200227.good_seq_ident.txt  # confirmed kinase PDB
            |------ check_kinase.200227.check_seq_ident.txt # ambigious kinase PDB, need check
            |------ check_kinase.200227.bad_seq_ident.txt   # unlikely kinase PDB
            |------ check_kinase.200227.no_seq_ident.txt    # confirmed non-kinase PDB
            |------ check_kinase.200227.non_kinase.txt      # checked non-kinase PDB
            |------ check_kinase.200227.checked_pdb.list    # all checked PDB
            |------ kinase_query_template.setup    # setup
```
- This example searches the RCSB PDB to first retreive the sequences that fit the kinase queries (all types of protein kinases) and compare them to canonical human kinome, which removes anything that falls show of 40% identity to any of the kinase sequences. It also compare to previously downloaded and checked PDB to avoid downloading the same thing or non-kinase again. Output are the multi-sequence alignment of retreived kinases to human kinome and their PDB structures (individual chains, save one unique kinase per PDB) superposed to the reference structure 1ATP_E.

##########################################################################
- **Example 2: Calculate structural parameters and Classify Kinase Conformation**
```
> examples/
    |--2_classify_kinase_conf
            |------ run_classifier.sh       # running this example
            |------ *.1atp.pdb              # superposed PDB
            |------ kinfo_pdb.list          # superposed PDB list
            |------ kinfo_pdb.fasta         # alignment of PDB seq
            |------ Template_parameter_file.in   # setup
            |
            |------ 1_result
                      |------ kinfo_pdb.example.csv  # kinase structural parameters
                      |------ kinfo_pdb.example.xlsx # kinase structural parameters
                      |------ kinfo_pdb.example.mssing.txt  # kinase with missing residues
                      |------ kinfo_pdb.example.SK_gb_kinfo_classify.csv # GB-based classification result
            |------ x_work
                      |------ kinfo_pdb.example.all_comb.fasta  # cleaned alignmnet of input
                      |------ _TEMP.temp_comb.fasta   # temporary alignment of input
```
- This examples runs both the structural parameter extraction from superposed PDB structures, and the Conformation Classification using the extracted parameters.

##########################################################################
- **Example 3: Run Conformation Classification only**
```
> examples/
    |--3_run_classifier_only
            |------ run_classifier_only.sh  # running this example
            |------ test.SK_*_kinfo_classify.csv  # result
            |------ output.200224.csv       # kinase structural parameters
```
- This example runs just the Conformation Classification portion and takes in only a pre-calculated set of structural parameters.

##########################################################################
- **Example 4: Get PDB structure info**
```
> examples/
    |--4_get_pdb_info
            |------ run_get_pdb_info.sh    # running this example
            |------ run_get_ligand_info.sh # running this example
            |------ test.list              # list of PDB to extract data
            |------ test-sample.csv.gz     # sample PDB result, csv
            |------ test-sample.xlsx       # sample PDBresult, xlsx
            |------ test-again.csv.gz      # sample ligand result, csv
            |------ test-again.xlsx        # sample ligand result, xlsx
```
- This example runs through all PDB structures (particular chain) and retreive the relevant data contained in PDB Header, e.g. publication date, ligands, mutations, modifications, resolutions, etc.; the other script will read from the file of Parsed PDB Header file to download the SMILES string of bound ligands and associate these ligands to the "pdb_id" they came from.

##########################################################################
- **Required packages**
```
  csh/tcsh        # Shell
  PyMOL           # 1.7.0+

  Blast+          # 2.2.6
    blastp

  muscle          # 3.8.31

  Python          # 3.6.8+
    biopython     # 1.72+
    pandas        # 0.24.2+
    numpy         # 1.16.2+
    sklearn       # 0.22.1  # may not be backward compatible
    tqdm          # 4.31.1+
    pathos        # 0.2.3+
    rdkit         # 2019.09.2
    requests      # 2.21.0
    time          #
    xmltodict     # 0.12.0
    tzlocal       # 2.0.0

Retired
#    r-randomforest # 4.6_14
#    rpy2          # 2.9.4   # known bug of not having "tzlocal"

```