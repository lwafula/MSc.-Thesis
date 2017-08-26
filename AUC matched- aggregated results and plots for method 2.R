
# aggregated AUC
# using only those actively-called treatments that appear in untransformed, cell- and well-level transformed data
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
#install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival") 
# source("http://bioconductor.org/biocLite.R") 
# biocLite(c("GO.db", "preprocessCore", "impute")) 

library(WGCNA) #https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/InstallationInstructions.html
library(knitr) #for kable
library("reshape2")

#run the run_mRMR function
setwd('E:/uhasselt/22yrsem/MThesis')

datawd<-'E:/uhasselt/22yrsem/MThesis/data'
datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/data/glog transformed datasets by lambda - Copy'
datawdmatrixpvals<-'E:/uhasselt/22yrsem/MThesis/report/results-matricesH2_pvalues etc'
datawdauc<-'E:/uhasselt/22yrsem/MThesis/report/dataView/aucmethod'

# initial dataset used for the untransformed data
p1601xxPPS_means_preprocessed<-'data/p1601xxPPS_means_preprocessed.Rdata'
load(p1601xxPPS_means_preprocessed); rm(p1601xxPPS_means_preprocessed)

#use only the active treatments
datp.orig<-  datp[datp[['Active_FractionOfRepl']]>=0.5,]
rm(datp)
lambv<-c(8.5,10,10.5,11)
for(i in 1:length(lambv)){
#cell-level transformed data
fileb1<-paste0("p1601xxPPS_means_preprocessedfeat_act_",lambv[i],'.RData')
fileb1<-file.path(datawdtranslambda,fileb1)
#lambda value 
lambdaval=lambv[i]; lambdaval

load(fileb1)
datpopt.featurelambda<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,]
rm(list=c('glogdatp','fileb1'))

# well level
fileb1<-paste0("aggregated_p1601xxPPS_means_preprocessedfeat_act_",lambv[i],'.RData')
fileb1<-file.path(datawdtranslambda,fileb1)
load(fileb1)
welldatpopt.featurelambda<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,]
rm(list=c('glogdatp','fileb1'))

# for all, just keep those treatments that were active in both untransformed and transformed [cell and well level] data sets
int.all<-Reduce(intersect, list(names(table(datp.orig$Treatment)),names(table(datpopt.featurelambda$Treatment)),
                                names(table(welldatpopt.featurelambda$Treatment)))) # 682 for lambda= 10.5

#datpOpt.feature.temp<-datpOpt.feature
datpOpt.feature<- datp.orig[datp.orig$Treatment %in% int.all,]
datpopt.featurelambda<- datpopt.featurelambda[datpopt.featurelambda$Treatment %in% int.all,]
welldatpopt.featurelambda<-welldatpopt.featurelambda[welldatpopt.featurelambda$Treatment %in% int.all,]

# AUC for untransformed
################### correction from Steffen
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
# optN <- FSres$optimalNFeatures # optimal number of features
mpya<-mpya[,c('value','nFeatures','run','gr')]

# cell level
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
mpya.1<-mpya.1[,c('value','nFeatures','run','gr')]
mpya<-rbind(mpya,mpya.1)

#well level 
nFeaturesConsidered <- length(attr(welldatpopt.featurelambda, "features_mRMR_optN"))
FSres.featurelambda <- evaluationAUC(
  datasetActive       = as.data.frame(welldatpopt.featurelambda),
  features           = attr(welldatpopt.featurelambda, "features_mRMR_optN"),
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

mpya.1$gr<-'2'
mpya.1<-mpya.1[,c('value','nFeatures','run','gr')]
mpya<-rbind(mpya,mpya.1)

## combine the plots for all AUC for all lambdas against untransformed
mpya$gr<-as.numeric(mpya$gr)
mpya$gr<-factor(mpya$gr,levels = c(0,1,2),labels = c('untransformed','cell level','well level'))
auc<-paste0('matched_cellwell/aggregated_auc',lambv[i],'.pdf')
auc<-file.path(datawdauc,auc)
pdf(auc)
plot(value~as.factor(gr),data=mpya,xlab='',
     ylab='AUC of replicates vs non-replicates across runs',ylim=c(min(mpya$value)-0.005,max(mpya$value)+0.005),main='Evolution of AUC assessing replicability Vs
     transformation parameter',pch=19)
abline(h=median(FSres$aucAll),col='red',lwd=2,lty=2)
legend('bottomleft',legend=c(expression(paste('median for untransformed'))),lty=c(2),col=c('red'),lwd=c(2),bty='n')
dev.off()
}












