
# step 1 -- retreive new kinases that haven't already been downloaded/examined
./1_find_new_kinase_pdb/run_find_new_kinase.sh

# step 2 -- check alignment and superposed kinases
jalview ./1_find_new_kinase_pdb/check_kinase.200227.fasta
pymol ./1_find_new_kinase_pdb/check_kinase.200227.super.pse

# step 3 -- calculate structural parameters and classify conformation
./2_run_kinfo_classifier/run_classifier.sh

  ** might need to run in 2 passes - 1st pass to assess all structure to find
     structural parameters that the script missed to identify, 2nd pass to 
     incorporate the missing info to finish the classification
     * use 3_make_correction_pymol.py

# step 4 -- get PDB info from retreived PDB
./4_get_pdb_info/run_get_pdb_info.sh
./4_get_pdb_info/run_get_ligand_info.sh
