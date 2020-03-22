
## this will initiate search of RCSB for new kinase xtal PDB
## it they are not already in the known kinsae xtal list
##  - all_downloaded_kinase_pdb.200227.list
## this will also skip any non-kinsae structure in the list
##  - check_kinase.200227.non_kinase.list


../../1_find_new_kinase_pdb.py \
  -r kinase_query.template.setup


## this will use the same setup script as before, but will not
## search the RCSB for every possible kinase xtal with queries.
## Instead it will only search for the listed Xtal PDB_ID

#../../1_find_new_kinase_pdb.py   \
#  -r kinase_query.template.setup \
#  -p known_kinase_set.list
