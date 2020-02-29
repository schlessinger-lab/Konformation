#!/usr/bin/env python3

import sys,os
import re
import bz2
import time
import pickle
import sklearn
import numpy as np
import pandas as pd

from argparse import ArgumentParser
from sklearn.impute import SimpleImputer

from x_konf_vars import DefaultVariables
from x_kinfo_SK_classify import KinfoClassify

np.seterr(invalid='ignore')

################################################
msg ='''
  > {0}
      -data   <file>      [ CSV file of re-calculated PDB structural parameters ]
      -mdldir <path>      [ Path to directory with ML models ]
      -out    <prefix>    [ Output prefix ]\n
  Optional:
      -use_sk <model>     [ Use SKLearn ML model: rf|svm|nn|kn|dt|gp|gb (def: rf) ]
      -use_r_rf           [ Use R::randomForest instead of SKLearn Classifier (def: None) ]\n
e.g.>  4_kinase_conf_ML_only.py
          -data data.csv  -use_sk gb
          -lib '/Users/xxx/scripts/Kinformation/z_database'
'''.format(sys.argv[0])
if len(sys.argv) == 1: sys.exit(msg)

##########################################################################

Selected_Cols = [ 'pdb_id','p1p1x','p2p2x','r3r3x','h_cgvc',
                  'ang_NHs','ang_CHs','dist_NH','dist_CH']

#########################################################################
def main():

  parm = DefaultVariables()
  args = UserInput()

  ## Check input SKL model option
  if args.use_sk is None:
    args.use_sk = 'rf'
    print('\033[31m INFO: No -use_sk supplied, fall back to default\033[0m: {0}'.format('rf'))
  elif args.use_sk not in ['rf','svm','nn','dt','kn','gb','gp']:
    sys.exit('\033[31m ERROR: -use_sk is invalid\033[0m: {0}'.format(args.use_sk))

###################

  ## Load pre-calculated kinase structural parameters
  df = pd.read_csv(args.data_file,sep=',',comment='#')[Selected_Cols]
  df.set_index('pdb_id', inplace=True)

  ## Impute any missing value in dataset before ML classication
  imputer = SimpleImputer(missing_values=np.nan, strategy='median')
  Imp_lst = imputer.fit_transform( df.to_numpy() )  # skip 1st column, pdb_id

  data_df = pd.DataFrame(Imp_lst, index=df.index, columns=df.columns)
  data_df.reset_index(inplace=True)

  ## use Kinformation Classifier to assign kinase conformation/confidence
  start = time.perf_counter()
  KinfoClassify(data_df, args.model_dir, args.outpref, parm,
                args.use_r_rf, args.use_sk)
  end   = time.perf_counter()
  print('\n## Total time to SKLearn \033[31m{0}\033[0m Classification: \033[31m{1:.3f}\033[0m ms for \033[34m{2}\033[0m frames'.format(
        args.use_sk, (end-start)*1000, len(data_df)))


##########################################################################

def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-data', dest='data_file', required=True,
                  help='Pre-calculated PDB structural parameters')
  p.add_argument('-out', dest='outpref', required=True,
                  help='Output prefix')
  p.add_argument('-mdldir', dest='model_dir', required=True,
                  help='Path to directory with ML models')

  p.add_argument('-use_sk', dest='use_sk', required=False,
                  help='Use SKLearn ML model: rf|svm|nn|kn|dt|gp|gb (def: rf)')
  p.add_argument('-use_r_rf', action='store_true',
                  help='Use R::randomForest instead of SKLearn RFClassifier (def: None)')

  return p.parse_args()

##########################################################################
if __name__ == '__main__':
  main()

##########################################################################
#
# Peter M.U. Ung @ MSSM/Yale
#
