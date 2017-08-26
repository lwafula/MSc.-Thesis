
#DATASET2
#  AUC boxplots for all lambdas
# using only those actively-called treatments that appear untransformed and respective transformed data sets
rm(list=ls())
# UPDATED on 15-08-2017 (LWafula)
library(MRMRFS)
library(hexbin)
library(ggplot2)
library(matlab)
library(plyr)
library(dplyr)
library(data.table)
library(flux)
library(hcs) 
library(WGCNA) #https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/InstallationInstructions.html
library(knitr) #for kable
library("reshape2")
#
setwd('E:/uhasselt/22yrsem/MThesis/DATASET2')

datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy'
datawdmatrixpvals<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/results-matricesH2_pvalues etc'
datawdauc<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/dataView/aucmethod'

# initial dataset used for the untransformed data
p1602xxPPS_means_preprocessed<-'E:/uhasselt/22yrsem/MThesis/data/p1602xxPPS_means_preprocessed.Rdata'
load(p1602xxPPS_means_preprocessed); rm(p1602xxPPS_means_preprocessed)

#use only the active treatments
datp.orig<-  datp[datp[['Active_FractionOfRepl']]>=0.5,]
rm(datp)
lambv<-c(0.1, seq(0.5,25,0.5))
lambdamatched<-c()
for(i in 1:length(lambv)){
  # cell-level transformed data
  fileb1<-paste0("p1602xxPPS_means_preprocessedfeat_act_",lambv[i],'.RData')
  fileb1<-file.path(datawdtranslambda,fileb1)
  #lambda value 
  lambdaval=lambv[i]; lambdaval
  load(fileb1)
  datpopt.featurelambda<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,]
  rm(list=c('glogdatp','fileb1'))
  
  # for all, just keep those treatments that were active in both untransformed and cell-transformed data sets
  int.all<-Reduce(intersect, list(names(table(datp.orig$Treatment)),names(table(datpopt.featurelambda$Treatment)))) # 642 for lambda= 0.5
  datpOpt.feature<- datp.orig[datp.orig$Treatment %in% int.all,]
  datpopt.featurelambda<- datpopt.featurelambda[datpopt.featurelambda$Treatment %in% int.all,]
  
  # AUC for untransformed
  columnReplicateID   <- "Treatment"  # the column that contains the replicate information; all rows with the same value are considered replicates
  cellcount_ft            <- 'CellCount_AllRetained'
  min_cellcount       <- 100
  nFeatures               <- 75   # number of features to select (initially, the optimal number which is usually less than nFeatures will be determined below)
  N_REPLICATES            <- 0    # > 0: randomly sample n replicates per treatment (to reduce influence of compounds with large number of replicates)
  AUCeval.m <- 10
  AUCeval.nPairsReplicates <- 0.2 
  AUCeval.nPairsNonReplicates <- 0.2
  AUCeval.nPairsReplicatesMax <- 10000
  AUCeval.nPairsNonReplicatesMax <- 10000
  AUCeval.distMethod <- "pearson"
  AUCeval.nThresholdsAUC <- 100
  
  nFeaturesConsidered <- length(attr(datpOpt.feature, "features_mRMR_optN"))
  FSres <- evaluationAUC(
    datasetActive       = as.data.frame(datpOpt.feature),
    features           = attr(datpOpt.feature, "features_mRMR_optN"),
    nFeaturesConsidered =nFeaturesConsidered,
    columnReplicateID   = columnReplicateID,
    m                   = AUCeval.m,
    nPairsReplicates            = AUCeval.nPairsReplicates,
    nPairsNonReplicates       = AUCeval.nPairsNonReplicates,
    nPairsReplicatesMax       = AUCeval.nPairsReplicatesMax,
    nPairsNonReplicatesMax = AUCeval.nPairsNonReplicatesMax,
    distMethod                          = AUCeval.distMethod,
    nThresholdsAUC                    = AUCeval.nThresholdsAUC
  )
  ################
  mpya<-melt(FSres)
  mpya<-mpya[complete.cases(mpya$nFeatures),]
  mpya$gr<-'0'
  mpya$lambda<-lambv[i]
  mpya<-mpya[,c('value','nFeatures','run','gr','lambda')]
  
  # cell transformed data
  nFeaturesConsidered <- length(attr(datpopt.featurelambda, "features_mRMR_optN"))
  FSres.featurelambda <- evaluationAUC(
    datasetActive       = as.data.frame(datpopt.featurelambda),
    features           = attr(datpopt.featurelambda, "features_mRMR_optN"),
    nFeaturesConsidered =nFeaturesConsidered,
    columnReplicateID   = columnReplicateID,
    m                   = AUCeval.m,
    nPairsReplicates            = AUCeval.nPairsReplicates,
    nPairsNonReplicates       = AUCeval.nPairsNonReplicates,
    nPairsReplicatesMax       = AUCeval.nPairsReplicatesMax,
    nPairsNonReplicatesMax = AUCeval.nPairsNonReplicatesMax,
    distMethod                          = AUCeval.distMethod,
    nThresholdsAUC                    = AUCeval.nThresholdsAUC
  )
  
  mpya.1<-melt(FSres.featurelambda)
  mpya.1<-mpya.1[complete.cases(mpya.1$nFeatures),]
  mpya.1$gr<-'1'
  mpya.1$lambda<-lambv[i]
  mpya.1<-mpya.1[,c('value','nFeatures','run','gr','lambda')]
  mpya<-rbind(mpya,mpya.1)
  
  # keep them all
  lambdamatched<-rbind(lambdamatched,mpya)
}

#                    save data
dsave<-file.path(datawdauc,paste0('matched_cell/data2allauc.RData'))
lambdamatched$gr<-as.numeric(lambdamatched$gr)
save(lambdamatched,file=dsave)

#                  box plots
#http://www.r-graph-gallery.com/9-ordered-boxplot/
# reorder the lambdas from the most sensible to the most resistant. (mixing low and high treatments for the calculations)
new_order <- with(lambdamatched, reorder(lambda , value, mean , na.rm=T))

# Then I make the boxplot, asking to use the 2 factors : variety (in the good order) AND treatment :
pdf<-file.path(datawdauc,paste0('matched_cell/D2_matched_auc_all_lambdas.pdf'))
pdf(pdf)
par(mar=c(4.5,4.5,3,1))
myplot=boxplot(value ~ gr*lambda , data=lambdamatched, boxwex=0.75 , ylab='AUC of replicates vs non-replicates across runs',
               main='',xlab=expression(paste(lambda,' ','value')),
               col=c("black" , "red") ,  xaxt="n",pch=19,cex=0.25,
               ylim=c(min(lambdamatched$value)-0.0005,max(lambdamatched$value)+0.0005))
title("Evolution of AUC assessing replicability Vs transformation parameter",cex.main = 1,font.main= 2, col.main= "black")

# To add the label of x axis
strsplit(myplot$names , '\\.')->lambdas
my_lambdas<-c()
for(i in 1:length(lambdas)){
  if(length(lambdas[[i]])==2){
    my_lambdas[i]<-lambdas[[i]][2]
  }
  if(length(lambdas[[i]])==3){
    my_lambdas[i]<-paste0(lambdas[[i]][2],'.',lambdas[[i]][3])
  }
}
my_lambdas=my_lambdas[seq(1 , length(my_lambdas), 2)]; my_lambdas<-sort(as.numeric(my_lambdas))
axis(1, at = seq(1.5 , 2*length(lambv) , 2), labels = my_lambdas , tick=TRUE , cex=0.3)
for(i in seq(0.55,2*length(lambv),4)){ abline(v=i,lty=1, col="black")}
# Add a legend
legend("bottomright", legend = c("untransformed", "transformed"), col=c("black" , "red"),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  ncol = 2, inset = c(0.1, 0.1))
dev.off()

ggplot(aes(y = value, x = lambda, fill = gr,group=lambda), data = lambdamatched) + geom_boxplot()
library(ggplot2)
ggplot(data = lambdamatched, aes(x = lambda, y = value,group=lambda)) + 
  geom_boxplot(aes(fill = gr), width = 0.35) + 
  theme_bw()+ scale_x_discrete(name =expression(paste(lambda,' ','value')),limits=c(0.1,seq(0.5,25,0.5)))+
  theme(axis.text.x = element_text(color="black", size=5, angle=45))
axis(1,at=seq(2,102,2),labels=my_lambdas,cex.axis=0.7)
















