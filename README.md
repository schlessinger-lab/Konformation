# Kinformation
Classification of Protein Kinase (Catalytic Domain) Conformations


__Citation 1:__ [\*Ung PMU, \*Rahman R, Schlessinger A. Redefining the Protein Kinase Conformational Space with Machine Learning. Cell Chemical Biology (2018) 25(7), 916-924.](https://doi.org/10.1016/j.chembiol.2018.05.002)

__Citation 2:__ [\*Rahman R, \*Ung PMU, Schlessinger A. KinaMetrix: a web resource to investigate kinase conformations and inhibitor space. 
Nucleic Acids Research (2019) 47(D1), D361–D366.](https://doi.org/10.1093/nar/gky916)

##########################################################################
##########################################################################
##########################################################################
- Updating the PDB structure library of Newly Published Protein Kinase Catalytic Domain
```
5_update_kinase_db.py
  -set             [ generate a default setup file using -r ]
  -r <setup name>  [ read setup file, perform RCSB PDB search for new kinases ]
  -p <list name>   [ skip RCSB search, download the PDB of known kinase listed ]

 e.g.> *.py -r setup.file
```
- The purpose of this script is to search the RCSB Protein Data Bank for all kinase-associated structures with a search Queries, then compare the new query results (PDB_ID) to an existing internal kinase PDB library (known_list) to find those that aren't in the library yet. Since RCSB search query is 'dirty' and includes anything that has association to kinase (kinase subdmain, or any protein that has assocation with the name "kinase"), these new entries are checked by their FASTA sequence (seq length > 220 + seq identity > 40% to human kinome) before downloading the kinase structures. The downloaded structures are then superposed onto 1ATP using a standard reference residue set.

- The new results are added to the existing internal kinase library list (known_list) for future reference. The script will not check and download structures with PDB and chain ID that have been checked before to avoid redundent download.

- As long as the (known_list) is kept current and the downloaded PDBs are placed into the correct folder, this script can be made into a repeating "cron" job to update the kinase structure library.

- I tried using MPI to speed up the download of FASTA and PDB but seems to run into server lockdown if multiple jobs are retrieving from the same IP. Could be a defense mechanism against DOS attack?? Using serial retrieve with a delay of 0.1s in between seem to work fine. It took ~ 10 mins to finish the checking and downloading of hundreds of new structures after 2.5 yrs of no update. But if the update is kept up every few months, the process should be much faster.


```
6_update_Kinfo_db.py

```

##########################################################################
- Required packages
```
csh/tcsh        # Shell
PyMOL           # 1.7.0+

Blast+          # 2.2.6
  blastp

t-coffee        # +
clustalo        # 1.2.4+
muscle          # 3.8.31

Python          # 3.6.8+
  biopython     # 1.72+
  pandas        # 0.24.2+
  numpy         # 1.16.2+
  tqdm          # 4.31.1+
  pathos        # 0.2.3+
  time          # 
  requests      # 2.21.0
```