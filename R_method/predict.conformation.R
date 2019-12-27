library(randomForest)
library(clusterSim)
library(dplyr)

## plug-in for classifying DFG/Chelix conformations in app.R
## load Pre-built RandomForest models
setwd('/Users/pmung/Dropbox (Schlessinger lab)/9_scripts/3_program/structures/4_Konformation/z_database')
#training_chelix_rf = get(load("R_rf_model_chelix.190527.rda")) # remove p1p1x,p2p2x,r3r3x data
#training_dfg_rf = get(load("R_rf_model_dfg.190527.rda")) # added r3r3x data

dat = read.table("stdy_kinase.param.171009.csv", header = T, sep ="," )

##get relevant columns 
dat_rel_col = dat[,c(1,4,6,7,8,9,10,12,13,14,15)]

##partition and format test/training data 
training = dat_rel_col[c(1,3:326),] ## row 1 is duplcate of row 2 
row.names(training) = training$pdb_id
training_class = as.factor(as.character(training$Group))
training_dfg = training$dfg_st
training[,c(1,2,10)] = NULL

test = dat_rel_col[327:nrow(dat_rel_col),]
row.names(test) = test$pdb_id
test[,c(1,2,10)] = NULL

##remove row with NA lose about 139 cases 
test_complete = test[complete.cases(test),]

##impute missing training data using c-helix class 
#training.impute = rfImpute(x = training, y = training.class ,   ntree = 1000)
#training.impute$training.class = NULL
##normalize data to 0
colMax <- function(data) sapply(data, max, na.rm = TRUE)

combined_data = as.data.frame(rbind(training, test_complete))
#combined_n = data.Normalization(combined_data[c(4,5,6,7)], type = "n5", normalization = "column" ); head(combined_n)
combined_n = combined_data
c_row = nrow(combined_data)

combined_part = combined_data[c(4,5,6,7)]; nrow(combined_part)

maxFunc = function(var, cb_mean) {  return(max(abs(var - cb_mean)))  }
cb_mean = apply(combined_part, 2, mean)
cb_max  = as.data.frame(Map(maxFunc, combined_part, cb_mean))
#cb_mean = as.data.frame(lapply(apply(combined_part, 2, mean), rep, nrow(combined_part)))
#cb_max  = as.data.frame( lapply( Map(maxFunc, combined_part, cb_mean), rep, c_row))
head(cb_mean); head(cb_max)


training_n = combined_n[1:325,]
test_n = combined_n[326:nrow(combined_n),]
head(test_n)

predict_conformation = function(desc) {
  
  desc_pred = desc[2,]
  row.names(desc_pred) = desc_pred$pdb_id
  desc_pred = desc_pred[,c(3,4,5,7,8,9,10,6)]  
  
  
#  combined_data_2 = combined_data = as.data.frame(rbind(combined_data, desc_pred))
#  combined_data_2 = as.data.frame(rbind(combined_data, desc_pred))

  ## Normalization of some columns using clusterSim function
  #combined_n = data.Normalization(combined_data_2, type = "n5", normalization = "column" )
  ## in case clusterSim is not available, use formula
#  combined_n = combined_data_2
  #combined_n = data.Normalization(combined_data[c(4,5,6,7)], type = "n5", normalization = "column" )

  ## normalization: if desc_pred < max of example set, use that, if > max, 
  normFunc = function(var, cb_mean, cb_max) { return((var - cb_mean)/cb_max) }
  temp = as.data.frame(Map(normFunc, desc_pred, cb_mean, cb_max)); head(temp)
  combined_n[c('ang_NHs','ang_CHs','dist_NH','dist_CH')] = temp
  
  desc_n = combined_n[nrow(combined_n),]
  test_pred_dfg = as.data.frame(
    predict( object = get(load("R_rf_model_dfg.rda")), 
             newdata = desc_n, type = "prob"))
  
  col = (which.max(test_pred_dfg[1,]))
  dfgstat = as.numeric(row.names(as.data.frame(col)))
  
  desc_n$dfg = as.factor(dfgstat)
  desc_n$dfg = as.numeric(as.character(desc_n$dfg))
  
  test_chelix_pred = as.data.frame(
    predict( object = get(load("R_rf_model_chelix.rda")),
             newdata = desc_n, type = "prob"))
  
  return(test_chelix_pred)

}