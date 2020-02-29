#!/usr/bin/env python3

import sys,os
import pandas as pd

import rpy2
import rpy2.robjects as ro
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()
R = ro.r

##########################################################################
## Run the R-implementation of Random Forest classification.
## Pre-trained models in .rda format can be loaded in
def R_RunRandomForest( data_df, lib_dir, models='' ):

  if not models:
    sys.exit('\033[31m  ERROR: No R randomForest model supplied\033[0m')
  else:
    R_dfg_model, R_chx_model = models

    if not os.path.isfile(lib_dir+R_dfg_model):
      sys.exit('\033[31m  ERROR: R randomForest model not found:\033[0m '+R_dfg_model)
    if not os.path.isfile(lib_dir+R_chx_model):
      sys.exit('\033[31m  ERROR: R randomForest model not found:\033[0m '+R_chx_model)

    Rfunc = R( """
function( dfg_mod, chx_mod ) {
  rf_model.dfg = load(dfg_mod)
  rf_model     = load(chx_mod)

  return( list(rf_model.dfg, rf_model) ) }  """ )

    models = Rfunc(lib_dir+R_dfg_model, lib_dir+R_chx_model)

#####################
  Rfunc = R( """ 
function( test, rf_model.dfg, rf_model ) {
  library(randomForest)

  ## rf-predict DFG status for test set 
  test.rf_pred.dfg = as.data.frame(
      predict( object = rf_model.dfg, newdata = test, type = "prob") )
#  head(test.rf_pred.dfg) ; str(test.rf_pred.dfg)

  ## assign rf-determined DFG conformation into original test set
  test.dfg = test 
  for (i in 1:nrow(test.rf_pred.dfg)) {
    col = which.max(test.rf_pred.dfg[i,])
    dfgstat = as.numeric(row.names(as.data.frame(col)))
    test.dfg[i, 9] = dfgstat 
  }

  ## change column name V9 to 'dfg_conf', and make factorize it
  colnames(test.dfg) = c(colnames(test), 'dfg_conf') #; head(test.dfg)
  test.dfg$dfg_conf = factor(test.dfg$dfg_conf)

  ## Create a new df, Add a column Row.names with PDB ID
  test.dfg.pred = test.dfg #; head(test.dfg.pred)
  test.dfg.pred$Row.names = row.names(test.dfg)
  test.dfg.pred[,c(1:8)] = NULL #; head(test.dfg.pred)

  ##predict probablities for classes for test set based on training 
  test.pred = as.data.frame(
      predict( object = rf_model, newdata = test.dfg, type = "prob") )
  test.pred.class = test.pred #; head(test.pred.class)

  ##assign class to rows based on class with max probablity 
  for ( i in 1:nrow(test.pred)) {
    col = which.max(test.pred[i,])
    test.pred.class[i,6] = as.character(row.names(as.data.frame(col)))
    test.pred.class[i,7] = as.numeric(test.pred[i,col])  
  }
  colnames(test.pred.class) = c(colnames(test.pred), "Class", "Probility")
  test.pred.class[,1:5] = NULL #; head(test.pred.class)

  ## combine all data, ignore those with NAs by associating using row.names
  predicted.class = merge(test.pred.class, test, by='row.names')
  head(predicted.class) #; nrow(predicted.class)
  predicted.class = merge(predicted.class, test.dfg.pred, by = 'Row.names')
  colnames(predicted.class)[1] = 'pdb_id'

  return( predicted.class )
} """ )

  R_rst   = Rfunc(pandas2ri.py2ri(data_df), models[0], models[1])
  R_rf_df = pandas2ri.ri2py( R_rst )
#  R_rf_df.set_index('pdb_id')
  return R_rf_df


########################################################################
########## Only used for first time to generate R RF model #############
########################################################################
## R-randomForest rfImpute function, takes categorized "responses" to
## subset the data and impute each subset with randomForest, hence all
## imputed data is slightly different and look more 'natural' than with
## SK Imputer, which has same imputed value for all missing data.
def R_Impute( df, Kinfo ):

  Rfunc = R( """ 
function( rdf, kinfo ) {
  library(randomForest)
  kinfo  = factor(kinfo)      # categorize $Group for rf
  rdf_impute = rfImpute(x = rdf, y = kinfo, ntree = 100)
  rdf_impute$kinfo = NULL

  return(rdf_impute)    # Remove categorized columns to avoid bug
  }  """ )

  df['kinfo'] = Kinfo

  # bug with rpy2, cannot handle string iput well
  raw = df.reset_index()
  col = ['pdb_id','Group','dfg_conf','kinfo']
  rdf = ro.conversion.py2ri( raw.drop(col, axis=1) )
  kin = ro.conversion.py2ri( raw.kinfo )

  impute = pandas2ri.ri2py( Rfunc( rdf, kin ) )

  ## Need to remove categorized columns due to rpy2 v2.9 bug. Should be
  ## fixed by v3.0 and won't need to remove those columns on R side 
  train_df = pd.concat( [df['Group'], impute.set_index(df.index)], axis=1 )
  return train_df


#################################################################
## Train Random Forest models with R
def R_TrainRandomForest( train_df, lib_dir ):

  Rfunc = R( """  
function( complete_set, lib_dir ) {
  library(randomForest)
#  library(clusterSim)    ## conda has no clusterSim

  Maximum   = function(var, cb_mean) { return(max(abs(var-cb_mean))) }
  Normalize = function(var, cb_mean, cb_max) { return((var-cb_mean)/cb_max) }

  ## normalize ang_ and dist_ data but leave vector data unchanged
  ## use data from both test and training sets to get wider range of input
#  complete_set  = as.data.frame(rbind(train_x, test_x))
  data_pre_norm = complete_set[c('ang_NHs','ang_CHs','dist_NH','dist_CH')]
  d_row = nrow(complete_set)

  ## create a repeating list of mean values for normalization vectorization
  cb_mean = as.data.frame( lapply( apply(data_pre_norm, 2 mean) ),rep, d_row )
  cb_max  = Map(Maximum, data_pre_norm, cb_mean)
  cb_max  = as.data.frame(lapply(cb_max, rep, d_row))
  temp    = as.data.frame(Map(Normalize, data_pre_norm, cb_mean, cb_max))
  complete_set[c('ang_NHs','ang_CHs','dist_NH','dist_CH')] = temp

  ## split the normalized dataset back into training and testing sets
  train = complete_set[1:nrow(train_x),]
  train$dfg_conf = factor(train$dfg_conf)
  train$Group    = factor(train$Group)

  ## train DFG predictor for: DFG-in, -out, -intermediate(other)
  ## include r3r3x and exclude h_cgvc + ang_ data, best OOB ~ 3.69%
  dfg_train_data = train['p1p1x','p2p2x','dist_NH','dist_CH','r3r3x')]
  dfg_train_conf = train$dfg_conf
  rf_model_dfg = randomForest(dfg_train_conf ~., data = dfg_train_data,
                              ntree = 1000)
#  rf_model_dfg$importance; rf_model_dfg
  save(rf_model_dfg, file = 'R_rf_model_dfg.rda')
#  write.table(rf_model_dfg$confusion, file = 'R_rf_model_dfg.confuse.txt', sep='\t')
#  write.table(rf_model_dfg$confusion, file = 'R_rf_model_dfg.importance.txt', sep='\t')

  ## train Chelix predictor, use DFG_conf data from DFG Random Forest model
  rf_train_data = train(c('h_cgvc','ang_NHs','ang_CHs','dist_NH','dist_CH','dfg_conf'))
  rf_train_conf = train$Group
  rf_model = randomForest(train$Group ~., data = rf_train_data, ntree = 1000)
#  rf_model$importance ; rf_model
  save(rf_model, file = 'R_rf_model_full.rda')
#  write.table(rf_model_full$confusion, file = 'R_rf_model_full.confuse.txt', sep='\t')
#  write.table(rf_model_full$confusion, file = 'R_rf_model_full.importance.txt', sep='\t')

  return( list(rf_model_dfg, rf_model) )
 }  """ )

  trained = Rfunc( pandas2ri.py2ri(train_df), lib_dir )
  return trained


#########################################################################
