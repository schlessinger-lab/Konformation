#!/usr/bin/env python3

import sys,os
import re
import bz2
import time
import pickle
import numpy as np
import pandas as pd

from collections import Counter
from sklearn.impute import SimpleImputer

from x_konf_vars import SKLearnDFGModelFiles
from x_konf_vars import SKLearnKinfoModelFiles
#from x_kinfo_R_classify import R_RunRandomForest  # rpy2 isn't well maintained

###################

Ref_Test_Cols = [ 'pdb_id','p1p1x','p2p2x','r3r3x','h_cgvc',
                  'ang_NHs','ang_CHs','dist_NH','dist_CH']

Ref_Final_Cols= [ 'pdb_id','Class',
                  'cidi_prob','cido_prob','codi_prob','codo_prob','wcd_prob',
                  'dfg_conf','dfg_prob','p1p1x','p2p2x','r3r3x','h_cgvc',
                  'ang_NHs','ang_CHs','dist_NH','dist_CH']

norm_cols = ['ang_NHs','ang_CHs','dist_NH','dist_CH']

dfg_train_cols = ['p1p1x','p2p2x','r3r3x','dist_NH','dist_CH']
full_train_cols= ['h_cgvc','ang_NHs','ang_CHs','dist_NH','dist_CH','dfg_conf']


##########################################################################
## Process the collected structural data by Normalization against a known
## kinase PDB derived Normalization parameter (mean and max), then run thru
## the RandomForest model

def KinfoClassify( data_df, lib_dir, outpref, parm, use_sk='rf' ):

  ## make sure the input dataframe has same columns as RF models
  if data_df.columns.isin(Ref_Test_Cols).sum() != len(Ref_Test_Cols):
    print('  \033[31mERROR: Column in Trajectory not matching required: \033[0m')
    print(data_df.columns)
    print(Ref_Test_Cols)
    sys.exit()

  ## Load factors for data normalization, then Normalize ang_/dist_ data
  norm_param = pd.read_csv(lib_dir+parm['kinfo_norm_param'],sep=',',comment='#')
  data_df[norm_cols] = Normalization(data_df[norm_cols], norm_param=norm_param)


#############
  ## use R-generated RandomForest model for classification
#  if use_r_rf:
#    print('\033[34m## Loading R RandomForest models...\033[0m')
#    models = [ parm['R_dfg_model'], parm['R_chx_model'] ]
#    result_df = R_RunRandomForest(data_df, lib_dir, models=models)
#    result_df.to_csv(outpref+'.R_rf_kinfo_classify.csv', sep=',')
#    print('\n\033[34m Write to:\033[0m {0}{1}'.format(outpref+'.R_rf_kinfo_classify.csv'))
#    return None

##############
  ## Use SK-generated RandomForest model for classification
  if use_sk:

    sk_dfg_models = SKLearnDFGModelFiles()
    sk_chx_models = SKLearnKinfoModelFiles()

    ## Load SK ML models
    print('\033[34m## Loading trained SK ML models... \033[31m{0}\033[0m'.format(use_sk))
    rfc_dfg  = pickle.load(bz2.open(lib_dir+sk_dfg_models[use_sk], 'rb'))
    rfc_full = pickle.load(bz2.open(lib_dir+sk_chx_models[use_sk], 'rb'))

    ## Run SK ML models
    result_df = SK_RunML(data_df, use_sk, models=[rfc_dfg, rfc_full])
#    print(Counter(result_df.Class))
    print('\n\033[34mConformation   Counts\033[0m')
    for conf, num in Counter(result_df.Class).most_common():
      print('\033[35m     {0}\t\033[31m{1}\033[0m'.format(conf, num))

    ## Print out results of classification along with probability
    result_df.to_csv(outpref+'.SK_{0}_kinfo_classify.csv'.format(use_sk), sep=',', index=None)
    print('\n\033[34m Write to:\033[0m {0}.SK_{1}_kinfo_classify.csv'.format(outpref, use_sk))
    return None



##########################################################################

##########################################################################
## Run Classifier, input data must have no NaN value or SKL will fail
## also provide probability for each classification type
def SK_RunML( df, ml_alg, models='' ):

  rfc_dfg, rfc_full = models
  ##### classify DFG conformation of trajectory frames #####
  start = time.perf_counter()
  data_dfg_pred = rfc_dfg.predict(df[dfg_train_cols])
  data_dfg_prob = rfc_dfg.predict_proba(df[dfg_train_cols])

  # append 'dfg_conf' and probability data to traj frame data
  df['dfg_conf'] = data_dfg_pred
  df['dfg_prob'] = np.max(data_dfg_prob, axis=1)
  print(' \033[34mSK_RF Classify DFG:\033[0m   {:.6f} s'.format((time.perf_counter()-start)))

  ##### classify Chelix/DFG conformation of traj frames #####
  start = time.perf_counter()
  traj_full_pred = rfc_full.predict(df[full_train_cols])
  traj_full_prob = rfc_full.predict_proba(df[full_train_cols])
  print(' \033[34mSK_RF Classify Kinfo:\033[0m {:.6f} s'.format((time.perf_counter()-start)))

  ## append 'Class' and probability to traj frame data 
  start = time.perf_counter()
  df['dfg_conf']    = state_dfg(pd.DataFrame(data_dfg_pred))
  df['Class']       = state_kinfo(pd.DataFrame(traj_full_pred))
  df['cidi_prob'] = traj_full_prob[:,0]
  df['cido_prob'] = traj_full_prob[:,1]
  df['codi_prob'] = traj_full_prob[:,2]
  df['codo_prob'] = traj_full_prob[:,3]
  df['wcd_prob']  = traj_full_prob[:,4]
  print(' \033[34mAppend Kinfo data:\033[0m    {:.6f} s'.format((time.perf_counter()-start)))

  return df[Ref_Final_Cols]


########################################################################
## Normalize ang_ and dist_ data with vectorization, same result as R's
## clusterSim data.Normalization(input, type='n5', normalization='column')
def Normalization( data, norm_param='' ):

  if not len(norm_param):
    cb_mean = data.to_numpy().mean(axis=0)
    cb_vars = data.to_numpy() - cb_mean
    cb_max  = np.max(np.abs(cb_vars), axis=0)
  else:
    cb_vars = data.to_numpy() - norm_param['mean'].to_numpy()
    cb_max  = norm_param['max'].to_numpy()

  return cb_vars/cb_max  # (var-mean)/max(abs(var-mean))


########################################################################
## set DFG conformation type based on DFG in/out with vectorized T/F
## Pandas based and use half with numpy vectorization to give ~ 10x speedup
## old half-vectorized ~ 1.2s for 3800 items, full-vectorized ~5ms, ~240x speedup
def state_dfg( state ):
  conf_di = (state == 0)
  conf_do = (state == 1)
  conf = pd.DataFrame({ '0': ['other']*len(state) })
  conf[conf_di[0].to_numpy() == True] = 'di' # any DI is '0'
  conf[conf_do[0].to_numpy() == True] = 'do' # any DO is '1'
  return conf['0'].to_numpy()

################
def state_kinfo( state ):
  conf_cidi = (state == 0)
  conf_cido = (state == 1)
  conf_codi = (state == 2)
  conf_codo = (state == 3)
  conf = pd.DataFrame({ '0': ['wcd']*len(state) })
  conf[conf_cidi[0].to_numpy() == True] = 'cidi'  # ci DI is '0'
  conf[conf_cido[0].to_numpy() == True] = 'cido'  # ci DO is '1'
  conf[conf_codi[0].to_numpy() == True] = 'codi'  # co DI is '2'
  conf[conf_codo[0].to_numpy() == True] = 'codo'  # co DO is '3'
  return conf['0'].to_numpy()

########################################################################
