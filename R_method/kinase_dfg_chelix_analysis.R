###if library is not availble use the command install.packages('<packagename>')
##Load Libraries 
library(randomForest)
library(plotly)
##read file 
## set row names to be pdb_id
raw_data = read.table("./Dropbox (Schlessinger lab)/9_scripts/3_program/structures/4_Konformation/z_database/stdy_kinase.param.170825.rrahman.csv", sep = ",", header = T)
row.names(raw_data) = raw_data$pdb_id
head(raw_data)
y = raw_data

## remove redundant pdb_id column since row names are pdb_id 
y$pdb_id = NULL
y$pdb_id.1 = NULL
y = y[,c(2,4,5,6,7,8,10,11,13)] ; head(y)
y.train = y[c(1,3:265),] ; head(y.train) ; nrow(y.train)

##get rows where there are assigned group class (get manually classified pdb_id ) 
#y.train = y[which(y$Group != ""),] ; head(y.train) ; nrow(y.train)

##set DFG status as last column 0= out, 1=in ,2 = other 
for (i in 1:nrow(y.train))
{
    if (y.train[i, ]$Group == "other")
    {
      y.train[i,10] = 2
      
    }
    else
    {
        g = substr(y.train[i, ]$Group, 3, 5 )
        if (g == "di") {
          y.train[i, 10] = 1
        }
        else {
          y.train[i, 10] = 0              
        }
    }
}
head(y.train) ; nrow(y.train)
colnames(y.train)[ncol(y.train)] = 'dfg_conf' ; head(y.train)

################################################################################

##factoring data tells R that the values are categorical data 
y.train$Group = factor(y.train$Group)
y.train$dfg_conf = factor(y.train$dfg_conf)

##set seed for randomforest 
set.seed(10)

## imput null/missing data using columns 2:9 to predict chelix classes 
y.train.rf.impute = rfImpute(
  x = y.train[,c(2:10)],  y = y.train$Group, ntree = 1000 )
head(y.train.rf.impute)

## rename 'Group' column name messed up by rfImpute 
colnames(y.train.rf.impute)[1] = 'Group'

## Use all columns to train a DFG classifier   
y.train.rf.dfg = randomForest( 
  dfg_conf ~., data = y.train.rf.impute[ ,c(2:ncol(y.train.rf.impute))], 
  ntree = 1000 )
head(y.train.rf.dfg) ; y.train.rf.dfg$confusion
save.image(y.train.rf.dfg, file = 'dfgmod.rda')

## Use all columns, train Chelix predictor 
y.train.rf = randomForest( Group ~.,data = y.train.rf.impute, ntree = 1000 )
y.train.rf$confusion
save.image(y.train.rf, file = 'chelixmod.rda')

###############################################################################

## subset data to test set with the remaining rows (rows with no manual classes)
y.test = y[ which(y$Group == "")[1]:nrow(y), ]
head(y.test) ; nrow(y.test) ; ncol(y.test)
## remove first column in test set (empty since no classes)
y.test[,c(1)]= NULL ; head(y.test) ; nrow(y.test) ; ncol(y.test)
## remove any rows with null data 
y.test = y.test[complete.cases(y.test), ] ; nrow(y.test)

## rf-predict DFG status for test set 
y.test.rf_pred.dfg = as.data.frame(
  predict( object = y.train.rf.dfg, newdata = y.test, type = "prob") )
head(y.test.rf_pred.dfg) ; str(y.test.rf_pred.dfg)

## assign rf-determined DFG conformation for test set
y.test.dfg = y.test 
for (i in 1:nrow(y.test.rf_pred.dfg)) {
    col = (which.max(y.test.rf_pred.dfg[i,]))
    dfgstat = as.numeric(row.names(as.data.frame(col)))
    y.test.dfg[i, 9] = dfgstat
}

## change column name V9 to 'dfg_conf', and make factorize it
colnames(y.test.dfg) = c(colnames(y.test), "dfg_conf") ; head(y.test.dfg)
y.test.dfg$dfg_conf = factor(y.test.dfg$dfg_conf)

## Create a new df, Add a column Row.names with PDB ID
y.test.dfg.pred = y.test.dfg; head(y.test.dfg.pred)
y.test.dfg.pred$Row.names = row.names(y.test.dfg)
y.test.dfg.pred[,c(1:8)] = NULL ; head(y.test.dfg.pred)

##predict probablities for classes for test set based from training 
y.test.pred = as.data.frame(
  predict( object = y.train.rf, newdata = y.test.dfg, type = "prob"))
y.test.pred.class = y.test.pred ; head(y.test.pred.class)

##assign class to rows based on class with max probablity 
for ( i in 1:nrow(y.test.pred)) {
  col = which.max(y.test.pred[i,])
  y.test.pred.class[i,6] = as.character(row.names(as.data.frame(col)))
  y.test.pred.class[i,7] = as.numeric(y.test.pred[i,col])
}
colnames(y.test.pred.class) = c(colnames(y.test.pred), "Class", "Probility")
head(y.test.pred.class)
y.test.pred.class[,1:5] = NULL ; head(y.test.pred.class)

## combine all data, ignore those with NAs by associating using row.names
y.predicted.class = merge(y.test.pred.class, y.test, by="row.names")
head(y.predicted.class) ; nrow(y.predicted.class)
head(y.test.dfg.pred)
y.predicted.class = merge(y.predicted.class, y.test.dfg.pred, by = "Row.names")
colnames(y.predicted.class)[1] = 'pdb_id'
head(y.predicted.class)

write.table(x = y.predicted.class, file = "../Desktop/5_schles/Kinase inhibitor data set/6.26.17.predicted.chelix.dfg.conformation.csv", 
            sep =",", quote = F, eol = "\n", row.names = F, col.names = T )
write.table(x= y.train, file = "../Desktop/5_schles/Kinase inhibitor data set/6.26.17.verified.chelix.dfg.conformation.csv", 
            sep = ",", quote = F, eol = "\n", row.names = T, col.names = T)

###############################################################################

##3D scatterplot code for rf-predicted test data: 
colors <- c('#4AC6B7', '#1972A4', '#965F8A', '#FF7070', '#C61951')
p = plot_ly(
        y.predicted.class, x = ~h_cgvc, y = ~p1p1x, z = ~p2p2x, color = ~Class, 
        colors = colors,  mode = 'markers', size = ~1-Probability, 
        marker = list(symbol = 'circle', sizemode = 'diameter',opacity = .9), 
        sizes = c(10, 40),
        text = ~paste('Id:', Row.names, '<br>Class:', Class, '<br>Probability:', 
                      Probability))

y.train$Row.names = row.names(y.train)
y.train$Probablity = rep( "v", nrow(y.train))

y3 = y.train[,c(12,1,11,2:10)]

## Add the rest of predicted data to the same plot
names(y3) = names(y.predicted.class)
#colnames(y3) =  colnames(y.predicted.class)
y.combined = rbind(y3, y.predicted.class)



###############################################################################
####subcluster stuff 

y.combined.split = split(y.combined, f = y.combined$Class)
y.combined.split.na.omit = na.omit(y.combined.split$other)
fit.kmeans = kmeans(y.combined.split.na.omit[,c(4,7,8)], centers = 3)
hist(fit.kmeans$cluster)
y.combined.split.na.omit$clust = as.factor(as.vector(fit.kmeans$cluster))

y =plot_ly(
  y.combined.split.na.omit, x = ~h_cgvc, y = ~p1p1x, z = ~p2p2x, 
  color = ~clust, 
  text = ~paste('Id:', Row.names, '<br>Class:', Class, '<br>Probability:', 
                Probability) )

y= plot_ly(
  y.combined, x = ~h_cgvc, y = ~p1p1x, z = ~p2p2x, color = ~Class, 
  colors = colors,  mode = 'markers', 
  marker = list(symbol = 'circle', sizemode = 'diameter',opacity = .9), 
                sizes = c(10, 40),
                text = ~paste('Id:', Row.names, '<br>Class:', 
                              Class, '<br>Probability:', Probability))

meta = read.table("../Desktop/5_schles/Kinase inhibitor data set/jun_2017/stdy_kinase.xtal.complete.170601.csv", header = T, sep =",")
meta.ligand =meta[which(meta$ligand!="NaN"),]
meta.apo = meta[which(meta$ligand=="NaN"),]
y.combined[1,12] = "1ATP" 
y.combined[2,12] = "1ATP" 
y.combined.2 = y.combined

require(plyr)
y.combined.2[c(3:nrow(y.combined)),12] = ldply(strsplit(y.combined[c(3:nrow(y.combined)),]$Row.names, "_"))$V1

y.combined.ligand = y.combined.2[which(y.combined.2$V12 %in% meta.ligand$pdb_id),]
y.combined.apo = y.combined.2[which(y.combined.2$V12 %in% meta.apo$pdb_id),]

y = plot_ly(y.combined.ligand, x = ~h_cgvc, y = ~p1p1x, z = ~p2p2x, 
            color = ~Class, colors = colors,  mode = 'markers', 
            marker = list(symbol = 'circle', sizemode = 'diameter',opacity = .9), 
            sizes = c(10, 40),
            text = ~paste('Id:', Row.names, '<br>Class:', 
                          Class, '<br>Probability:', Probability))

y.combined.ligand.cidi = y.combined.ligand[y.combined.ligand$Class == "cidi",]
y.combined.ligand.cido = y.combined.ligand[y.combined.ligand$Class == "cido",]
y.combined.ligand.codo = y.combined.ligand[y.combined.ligand$Class == "codo",]
y.combined.ligand.codi = y.combined.ligand[y.combined.ligand$Class == "codi",]
y.combined.ligand.other = y.combined.ligand[y.combined.ligand$Class == "other",]
y.combined.ligand.cidi[which(y.combined.ligand.cidi$p2p2x > .75 & y.combined.ligand.cidi$p1p1x > .6),13] = 1 
y.combined.ligand.cidi[which(y.combined.ligand.cidi$p2p2x > .75 & y.combined.ligand.cidi$p1p1x <= .6),13] = 2 
y.combined.ligand.cidi[which(y.combined.ligand.cidi$p2p2x <= .75 & y.combined.ligand.cidi$p1p1x > .6),13] = 3 
y.combined.ligand.cidi[which(y.combined.ligand.cidi$p2p2x <= .75 & y.combined.ligand.cidi$p1p1x <= .6),13] = 4

y.combined.ligand.cidi$cluster_id =  factor(y.combined.ligand.cidi$V13)
y.combined.ligand.cidi$V13 = NULL

y.combined.ligand.cidi.split = split(y.combined.ligand.cidi, f = y.combined.ligand.cidi$cluster_id)
cidi.1.min = as.numeric(as.vector(apply(y.combined.ligand.cidi.split$`1`,2, min)[c(4,7,8)]))
cidi.1.max = as.numeric(as.vector(apply(y.combined.ligand.cidi.split$`1`,2, max)[c(4,7,8)]))
cidi.2.min = as.numeric(as.vector(apply(y.combined.ligand.cidi.split$`2`,2, min)[c(4,7,8)]))
cidi.2.max = as.numeric(as.vector(apply(y.combined.ligand.cidi.split$`2`,2, max)[c(4,7,8)]))
cidi.3.min = as.numeric(as.vector(apply(y.combined.ligand.cidi.split$`3`,2, min)[c(4,7,8)]))
cidi.3.max = as.numeric(as.vector(apply(y.combined.ligand.cidi.split$`3`,2, max)[c(4,7,8)]))
cidi.4.min = as.numeric(as.vector(apply(y.combined.ligand.cidi.split$`4`,2, min)[c(4,7,8)]))
cidi.4.max = as.numeric(as.vector(apply(y.combined.ligand.cidi.split$`4`,2, max)[c(4,7,8)]))
cidi.1.mid = (cidi.1.max - cidi.1.min) / 2
cidi.2.mid = (cidi.2.max - cidi.2.min) / 2
cidi.3.mid = (cidi.3.max - cidi.3.min) / 2
cidi.4.mid = (cidi.4.max - cidi.4.min) / 2
distance  = function(x, mid) {
    d = ( (as.numeric(x[4])- mid[1])^2 + (as.numeric(x[7])-mid[2])^2 + 
          (as.numeric(x[8])-mid[3])^2 )^.5
    print(d)
}

y.combined.ligand.cidi.split$`1`$distance_from_cluster_center = as.numeric(as.vector(apply(y.combined.ligand.cidi.split$`1`,1 , FUN = function(x) distance(x, cidi.1.mid) ) ))
y.combined.ligand.cidi.split$`2`$distance_from_cluster_center = as.numeric(as.vector(apply(y.combined.ligand.cidi.split$`2`,1 , FUN = function(x) distance(x, cidi.2.mid) ) ))
y.combined.ligand.cidi.split$`3`$distance_from_cluster_center = as.numeric(as.vector(apply(y.combined.ligand.cidi.split$`3`,1 , FUN = function(x) distance(x, cidi.3.mid) ) ))
y.combined.ligand.cidi.split$`4`$distance_from_cluster_center = as.numeric(as.vector(apply(y.combined.ligand.cidi.split$`4`,1 , FUN = function(x) distance(x, cidi.4.mid) ) ))
y.combined.ligand.cidi = rbind(y.combined.ligand.cidi.split$`1`,y.combined.ligand.cidi.split$`2`,y.combined.ligand.cidi.split$`3`,y.combined.ligand.cidi.split$`4`)


y.combined.ligand.codi[which(y.combined.ligand.codi$p2p2x > .6 & y.combined.ligand.codi$p1p1x > .5),13] = 1 
y.combined.ligand.codi[which(y.combined.ligand.codi$p2p2x > .6 & y.combined.ligand.codi$p1p1x <= .5),13] = 2 
y.combined.ligand.codi[which(y.combined.ligand.codi$p2p2x <= .6 & y.combined.ligand.codi$p1p1x > .5),13] = 3 
y.combined.ligand.codi[which(y.combined.ligand.codi$p2p2x <= .6 & y.combined.ligand.codi$p1p1x <= .5),13] = 4
y.combined.ligand.codi$V13 = factor(y.combined.ligand.codi$V13 ) 
y.combined.ligand.codi$cluster_id =  factor(y.combined.ligand.codi$V13)
y.combined.ligand.codi$V13 = NULL

y.combined.ligand.codi.split = split(y.combined.ligand.codi, f = y.combined.ligand.codi$cluster_id)
codi.1.min = as.numeric(as.vector(apply(y.combined.ligand.codi.split$`1`,2, min)[c(4,7,8)]))
codi.1.max = as.numeric(as.vector(apply(y.combined.ligand.codi.split$`1`,2, max)[c(4,7,8)]))
codi.2.min = as.numeric(as.vector(apply(y.combined.ligand.codi.split$`2`,2, min)[c(4,7,8)]))
codi.2.max = as.numeric(as.vector(apply(y.combined.ligand.codi.split$`2`,2, max)[c(4,7,8)]))
codi.3.min = as.numeric(as.vector(apply(y.combined.ligand.codi.split$`3`,2, min)[c(4,7,8)]))
codi.3.max = as.numeric(as.vector(apply(y.combined.ligand.codi.split$`3`,2, max)[c(4,7,8)]))
codi.4.min = as.numeric(as.vector(apply(y.combined.ligand.codi.split$`4`,2, min)[c(4,7,8)]))
codi.4.max = as.numeric(as.vector(apply(y.combined.ligand.codi.split$`4`,2, max)[c(4,7,8)]))
codi.1.mid = (codi.1.max - codi.1.min) / 2
codi.2.mid = (codi.2.max - codi.2.min) / 2
codi.3.mid = (codi.3.max - codi.3.min) / 2
codi.4.mid = (codi.4.max - codi.4.min) / 2

y.combined.ligand.codi.split$`1`$distance_from_cluster_center = as.numeric(as.vector(apply(y.combined.ligand.codi.split$`1`,1 , FUN = function(x) distance(x, codi.1.mid) ) ))
y.combined.ligand.codi.split$`2`$distance_from_cluster_center = as.numeric(as.vector(apply(y.combined.ligand.codi.split$`2`,1 , FUN = function(x) distance(x, codi.2.mid) ) ))
y.combined.ligand.codi.split$`3`$distance_from_cluster_center = as.numeric(as.vector(apply(y.combined.ligand.codi.split$`3`,1 , FUN = function(x) distance(x, codi.3.mid) ) ))
y.combined.ligand.codi.split$`4`$distance_from_cluster_center = as.numeric(as.vector(apply(y.combined.ligand.codi.split$`4`,1 , FUN = function(x) distance(x, codi.4.mid) ) ))
y.combined.ligand.codi = rbind(y.combined.ligand.codi.split$`1`,y.combined.ligand.codi.split$`2`,y.combined.ligand.codi.split$`3`,y.combined.ligand.codi.split$`4`)


y.combined.ligand.cido$cluster_id = rep(1,nrow(y.combined.ligand.cido))
cido.1.min = as.numeric(as.vector(apply(y.combined.ligand.cido,2, min)[c(4,7,8)]))
cido.1.max = as.numeric(as.vector(apply(y.combined.ligand.cido,2, max)[c(4,7,8)]))
cido.1.mid = (cido.1.max - cido.1.min) / 2
y.combined.ligand.cido$distance_from_cluster_center = 
  as.numeric(as.vector(apply(y.combined.ligand.cido,1 , FUN = function(x) distance(x, cido.1.mid) ) ))

y.combined.ligand.codo$cluster_id = rep(1,nrow(y.combined.ligand.codo))
codo.1.min = as.numeric(as.vector(apply(y.combined.ligand.codo,2, min)[c(4,7,8)]))
codo.1.max = as.numeric(as.vector(apply(y.combined.ligand.codo,2, max)[c(4,7,8)]))
codo.1.mid = (codo.1.max - codo.1.min) / 2
y.combined.ligand.codo$distance_from_cluster_center = 
  as.numeric(as.vector(apply(y.combined.ligand.codo,1 , FUN = function(x) distance(x, codo.1.mid) ) ))

y.combined.ligand.other$cluster_id = rep(1,nrow(y.combined.ligand.other))
other.1.min = as.numeric(as.vector(apply(na.omit(y.combined.ligand.other),2, min)[c(4,7,8)]))
other.1.max = as.numeric(as.vector(apply(na.omit(y.combined.ligand.other),2, max)[c(4,7,8)]))
other.1.mid = (other.1.max - other.1.min) / 2
y.combined.ligand.other$distance_from_cluster_center = 
  as.numeric(as.vector(apply(y.combined.ligand.other,1 , FUN = function(x) distance(x, other.1.mid) ) ))




y.combined.ligand = rbind(y.combined.ligand.cidi,y.combined.ligand.codo,y.combined.ligand.cido,y.combined.ligand.codi, y.combined.ligand.other)



write.table(y.combined.ligand, file = "../Dropbox (Schlessinger lab)/Schlessinger lab Team Folder/shared/1_kinase_family/ligand.bound.kinases.subcluster.distanace.6.28.17.csv", sep =",", quote = F, row.names = F , eol =  "\n")

plot_ly(y.combined.ligand.codi, x = ~h_cgvc, y = ~p1p1x, z = ~p2p2x, color = ~V13, 
        colors = colors,  mode = 'markers', 
        marker = list(symbol = 'circle', sizemode = 'diameter',opacity = .9), 
        sizes = c(10, 40),
        text = ~paste('Id:', Row.names, '<br>Class:', Class, '<br>Probability:', 
                      Probability))

kinase_param = read.table(
  "../Dropbox (Schlessinger lab)/Schlessinger lab Team Folder/shared/1_kinase_family/1_manual_classes/stdy_kinase.param.170626.csv",
  sep = "," , header = T )

kinase_param_fil = kinase_param[,c(1,15)]
colnames(kinase_param_fil) = c("Row.names","dfg_st")
y.combined.ligand.dfgst = merge(y.combined.ligand, kinase_param_fil, by = "Row.names")
y.combined.ligand.dfgst$dfg_st = factor(y.combined.ligand.dfgst$dfg_st)
y.combined.ligand.dfgst.sp = split(y.combined.ligand.dfgst, f= y.combined.ligand.dfgst$Class)
plot_ly(y.combined.ligand.dfgst,x= ~h_cgvc, y=~p1p1x ,z = ~p2p2x, color = ~dfg_st,
        text = ~paste('Id:', Row.names, '<br>Class:', Class, '<br>Probability:', 
                      Probability))


#klifs = read.table("../Desktop/5_schles/Kinase inhibitor data set/KLIF_classification_all_pdb_unique.csv", sep = ",")
#y.predictedy.predicted.class[which(y.predicted.class$Class != "other"), ] 


#api_create(p, filename = "predicted-coformations")
#heatmap.2(data.matrix(y.test.pred), trace = "none")
