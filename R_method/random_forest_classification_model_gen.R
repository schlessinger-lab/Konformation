library(tidyverse)
library(randomForest)
library(plotly)
library(clusterSim)
library(dplyr)

set.seed(10)
setwd('/Users/pmung/Dropbox (Schlessinger lab)/9_scripts/3_program/structures/4_Konformation')
dat = read.table("z_database/stdy_kinase.param.171009.csv",header =T, sep ="," )
head(dat)
##get relevant columns 
dat_rel_col = dat[,c(1,4,6,7,8,9,10,12,13,14,15)]
head(dat_rel_col)

##partition and format test/training data 
training = dat_rel_col[c(1,3:326),] ## row 1 is duplcate of row 2 
row.names(training) = training$pdb_id; head(training)
training_class = as.factor(as.character(training$Group))
training_dfg   = as.factor(training$dfg_st)
training[,c(1,2,10)] = NULL  # raw data
head(training)

test = dat_rel_col[327:nrow(dat_rel_col),]
row.names(test) = test$pdb_id; head(test)
test[,c(1,2,10)] = NULL; head(test)
##remove row with NA lose about 139 cases 
test_complete = test[complete.cases(test),]
head(test_complete)

##impute missing training data using c-helix class 
#training.impute = rfImpute(x = training, y = training.class ,   ntree = 1000)
#training.impute$training.class = NULL
##normalize data to 0
combined_data = as.data.frame(rbind(training, test_complete))
#combined_n = data.Normalization(combined_data[c(4,5,6,7)], type = "n5", normalization = "column" ); head(combined_n)

combined_n = combined_data
combined_part = combined_data[c(4,5,6,7)]; nrow(combined_part)
#cb_mean = apply(combined_part, 2, mean); cb_mean
#cb_max  = as.data.frame(Map(maxFunc, combined_part, cb_mean)); cb_max
cb_mean = as.data.frame(lapply(apply(combined_part, 2, mean), rep, nrow(combined_part)))
head(cb_mean)
maxFunc = function(var, cb_mean) {  return(max(abs(var - cb_mean)))  }
cb_max  = (Map(maxFunc, combined_part, cb_mean)); head(cb_max)
cb_max  = as.data.frame(lapply(cb_max, rep, nrow(combined_part)))
head(cb_mean); head(cb_max)

normFunc = function(var, cb_mean, cb_max) { return( (var - cb_mean) / cb_max ) }
temp = as.data.frame(Map(normFunc, combined_part, cb_mean, cb_max)); head(temp)
combined_n[c('ang_NHs','ang_CHs','dist_NH','dist_CH')] = temp
head(combined_n)

training_n = combined_n[1:325,]
test_n = combined_n[326:nrow(combined_n),]
head(test_n)


##fix dfg stat values
training_dfg_2 = c()

for (i in 1:length(training_class)) {
  if (training_class[i] == "other")  {
    training_dfg_2[i] = 3  }
  else  {
    g = substr(training_class[i], 3, 5 )
    if (g == "di") {
      training_dfg_2[i] = 1 }
    else {
      training_dfg_2[i] = 0 }  }  
}
training_dfg_2 = as.factor(training_dfg_2); head(training_dfg_2)

###train DFG classifier, include r3r3x, exclude ang_ data, best OOB ~ 3.69-4.0%
head(training_n[c(1,2,6:8)])
training_dfg_rf = randomForest( training_dfg_2 ~., data=training_n[c(1,2,6:8)], 
                                ntree = 1000); training_dfg_rf
training_dfg_rf$confusion
training_dfg_rf$importance
training_dfg_rf

##train Chelix classifer, remove p1p1x, p2p2x, r3r3x, best OOB ~2.46%
training_n_dfg = training_n
training_n_dfg$dfg = training_dfg_2 ##add dfg data to data 
head(training_n_dfg)
head(training_n_dfg[,c(3:7,9)])
training_n_dfg$dfg = as.numeric(as.character(training_n_dfg$dfg )); head(training_n_dfg)
training_chelix_rf = randomForest( training_class ~., 
                                   data=training_n_dfg[c(3:7,ncol(training_n_dfg))], 
                                   ntree = 1000); training_chelix_rf
training_chelix_rf$confusion
training_chelix_rf$importance

save(training_dfg_rf, file = "dfgmod.rda")
save(training_chelix_rf, file = "chelixmod.rda")
#xxx = get(load('chelixmod.rda'))   # need to use get() to complete retreive
#yyy = get(load('dfgmod.rda'))


##predict test DFG classes 
test_pred_dfg = as.data.frame(
                  predict(object = training_dfg_rf, newdata = test_n, 
                          type = "prob")) ; (test_pred_dfg)
test_pred_dfg_2 = test_n; head(test_pred_dfg_2)

for (i in 1:nrow(test_pred_dfg)) {
  col = which.max(test_pred_dfg[i,]); col
  dfgstat = as.numeric(row.names(as.data.frame(col)))
  test_pred_dfg_2[i, 9] = dfgstat  }
head(test_pred_dfg_2)

test_dfg =  test_pred_dfg_2$V9
test_pred_dfg_2$dfg = test_dfg
test_pred_dfg_2$V9 = NULL
head(test_pred_dfg_2)
test_chelix_pred = as.data.frame( 
  predict( object = training_chelix_rf, newdata = test_pred_dfg_2, type = "prob"))
test_chelix_pred
test_chelix_predictions = test_pred_dfg_2       # new dataset

##assign class to rows based on class with max probablity 
for ( i in 1:nrow(test_chelix_pred)) {
  col = which.max(test_chelix_pred[i,])
  test_chelix_predictions[i,10] =  as.character(row.names(as.data.frame(col)))
  test_chelix_predictions[i,11] = as.numeric(test_chelix_pred[i,col])
}
test_chelix_predictions$Class = test_chelix_predictions$V10
test_chelix_predictions$Probility = test_chelix_predictions$V11
test_chelix_predictions[,c(9,10,11)] = NULL
head(test_chelix_predictions)

write.table(x = test_chelix_predictions, 
  file = "../../2_predicted_classes/8.29.17.predicted.chelix.dfg.conformation.csv",
  sep =",", quote = F, eol = "\n", row.names = T, col.names = T )


head(test_chelix_predictions)
head(training)
training_n$dfg = training_dfg_2; head(training_n)
training_n$Class = training_class; head(training_n)
training_n$Probility = rep(1, nrow(training)) ; head(training_n)

full_data = as.data.frame(rbind(training_n, test_chelix_predictions))
head(full_data)

source = c(rep("training",nrow(training)),rep("test",nrow(test_chelix_predictions)))
head(source)

full_data$source = source; head(full_data)

write.table(x = full_data, 
  file = "../../2_predicted_classes/8.29.17.predicted+verified.chelix.dfg.conformation.csv",
  sep =",", quote = F, eol = "\n", row.names = T, col.names = T )


