#!/usr/bin/env python3

import sys
import pandas as pd


##########################################################################
def BuildDataSet():
  return { 
           'pdb_id':None,  'pdb':None,
           'h_phi':None,   'h_psi':None,   'h_axvc':None,
           'h_norm':None,  'h_cgvc':None,  'h_scvc':None,
           'h_ca':None,    'h_cg':None,    'h_curv':None,
           'n_ca':None,    'c_ca':None,
           'n_phi':None,   'n_psi':None,
           'c_phi':None,   'c_psi':None,
           'dist_NH':None, 'dist_NC':None, 'dist_CH':None,
           'ang_NHs':None, 'ang_CHs':None,
           'r_medi':None,  'r_std':None,
           'g_resi':None,  'g_coord':None,
           'p1p1x':None,   'p2p2x':None,   'dfg_st':None,
           'v1v1x':None,   'v2v2x':None,   'dfg_vs':None,
           # d0_x/d1_x = p1/p2 without cross Ca-to-Ca vec
           'r3r3x':None,   'd0_x':None,    'd1_x':None,	
           # h_sc_x = cg_vec.sc_vec	h_cg_x = cg_vec.cb_vec
           'h_sc_x':None,  'h_cg_x':None,
           'r_curv':None,  'r_curv_m2':None,
           'c_curv':None,  'c_curv_m2':None,
           'lig_name':None,'lig_id':None,  
           'lig_ha':None,  'lig_multi':None,
           'pock_out':None,'pock_hinge':None,
           'pock_sugar':None, 'pock_up':None,
           'pock_I5':None, 'pock_I5in':None,
           'pock_I5dp':None,'pock_dfg':None,
           'pock_bp':None, 'lig_type':None,
         }


################################################
def Data2Pandas(Data, output):

  df = pd.DataFrame.from_dict(Data, orient='index')
  df.index.name = 'pdb_id'

  # Columns of data to print out 
  columns = ['pdb_id',
             'p1p1x','p2p2x','h_cgvc','h_scvc',
             'ang_NHs','ang_CHs','dist_NH','dist_CH','dfg_st','r3r3x',
             'v1v1x','v2v2x','g_resi',
             'lig_name','lig_id','lig_ha','lig_multi','lig_type',
             'pock_out','pock_hinge','pock_sugar','pock_up','pock_I5',
             'pock_I5in','pock_I5dp','pock_dfg','pock_bp']

  reorder = df[columns]

  # Output to CSV
  reorder.to_csv(output+'.csv', sep=',', na_rep='NA', encoding='utf-8',
                 float_format='%.5f', header=True, )

  # Output to Excel
  reorder.to_excel(output+'.xlsx', header=True , na_rep='NA')


##########################################################################
# 
#   Peter Man-Un Un @ MSSM
#
#   v1.0    17.03.11
#   v2.0    17.03.17    add more parameters
#   v2.1    17.04.28    add ax-cg normal vector, h-ax dot product, cg vector
#   v3.0    17.05.09	add more parameters
#   v4.0    17.08.18	add ligand-types, chean up R/C-spines
#   v5.0    17.10.07	add non-crossed DFG vectors
#   v6.0    17.12.21    only output most essential parameters
#
##########################################################################
