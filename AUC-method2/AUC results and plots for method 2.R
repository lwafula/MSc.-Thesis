
#DATASET 2
rm(list=ls())
# LWafula (25-06-2017) MRMRFS script practise
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

#run the run_mRMR function
setwd('E:/uhasselt/22yrsem/MThesis/DATASET2')

datawd<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data'
datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy/for H2'
datawdmatrixpvals<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/results-matricesH2_pvalues etc'
datawdauc<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/dataView/aucmethod'

# initial dataset used for the untransformed data
p1602xxPPS_means_preprocessed<-'E:/uhasselt/22yrsem/MThesis/data/p1602xxPPS_means_preprocessed.Rdata'
load(p1602xxPPS_means_preprocessed); rm(p1602xxPPS_means_preprocessed)

#use only the active treatments
datp<- datp[datp[['Active_FractionOfRepl']]>=0.5,]

#
################### AUC for the untransformed data
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

nFeaturesConsidered <- length(attr(datp, "features_mRMR_optN"))
FSres <- evaluationAUC(
  datasetActive       = as.data.frame(datp),
  features           = attr(datp, "features_mRMR_optN"),
  nFeaturesConsidered =nFeaturesConsidered,
  columnReplicateID   = columnReplicateID,
  m                   = AUCeval.m,
  nPairsReplicates            = AUCeval.nPairsReplicates,
  nPairsNonReplicates       = AUCeval.nPairsNonReplicates,
  nPairsReplicatesMax       = AUCeval.nPairsReplicatesMax,
  nPairsNonReplicatesMax = AUCeval.nPairsNonReplicatesMax,
  distMethod                          = AUCeval.distMethod,
  nThresholdsAUC                    = AUCeval.nThresholdsAUC,
  returnDetailsAUC=TRUE
)
untrans<-'D2_untransformed.pdf'
untrans<-file.path(datawdauc,untrans)
pdf(untrans)
plot(x = FSres)
dev.off()

#distributions for the untransformed cases (replicates vs non-replicates)
#replicates vs non-replicates
for (j in 1:10){
  plotname<-paste0('D2_replicatesvsnonreplicates_run_',j,'_for_untransformed','.pdf')
  plotname<-file.path(datawdauc,plotname)
  pdf(plotname)
  plot(x = FSres, 
       type = "densityDistances",
       run = j,
       feat = as.character(
         FSres$optimalNFeatures
       )
  )
  dev.off()
}

################
mpya<-melt(FSres)
mpya<-mpya[complete.cases(mpya$nFeatures),]
mpya$gr<-'0'; mpya<-mpya[,c('gr','run','nFeatures','value')]
# optN <- FSres$optimalNFeatures # optimal number of features

######### 
# here we do not limit the number of optimal features, just working with active treatments
lambv<-c()
file.names <- dir(datawdtranslambda, pattern ="p1602xxPPS_means_preprocessedfeat_act_")
for (i in 1:length(file.names)){
  fileb1<-file.path(datawdtranslambda,file.names[i])
  lambv[i]<-gsub('p1602xxPPS_means_preprocessedfeat_act_','',file.names[i]); lambv[i]<-gsub('.RData','',lambv[i]); lambv[i]
  
  load(fileb1)
  datpopt.featurelambda<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,]
  rm(list=c('glogdatp','fileb1'))
  
 #auc comparison
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

mpya.1$gr<-lambv[i]; mpya.1<-mpya.1[,c('gr','run','nFeatures','value')]

mpya<-rbind(mpya,mpya.1)

#####
#### run specific replicate vs non-replicates distribution plot
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
  nThresholdsAUC                    = AUCeval.nThresholdsAUC,
  returnDetailsAUC=TRUE
)

#usual AUC plot
plotname<-paste0('D2_auc_for_lambda',lambv[i],'.pdf')
plotname<-file.path(datawdauc,plotname)
pdf(plotname)
plot(x = FSres.featurelambda, type = "aucVsFeatures")
dev.off()

#replicates vs non-replicates
for (j in 1:10){
  plotname<-paste0('D2_replicatesvsnonreplicates_run_',j,'_for_lambda',lambv[i],'.pdf')
  plotname<-file.path(datawdauc,plotname)
  pdf(plotname)
  plot(x = FSres.featurelambda, 
       type = "densityDistances",
       run = j,
       feat = as.character(
         FSres.featurelambda$optimalNFeatures
       )
  )
  dev.off()
}
}

## combine the plots for all AUC for all lambdas against untransformed
mpya$gr<-as.numeric(mpya$gr)
auc<-'D2_auc.pdf'
auc<-file.path(datawdauc,auc)
pdf(auc)
plot(value~as.factor(gr),data=mpya,xlab=expression(paste(lambda,' ','value')),
ylab='AUC of replicates vs non-replicates across runs',ylim=c(0.95,0.965),main='Evolution of AUC assessing replicability Vs
transformation parameter',pch=19)
abline(h=median(FSres$aucAll),col='red',lwd=2,lty=2)
legend('bottomleft',
       legend=c(expression(paste(lambda,' ', '= 0 (untransformed case)')),expression(paste('median for',' ',lambda,' ', '= 0'))),
       pch=c(19,NA),lty=c(NA,2),col=c(NA,'red'),lwd=c(NA,2),bty='n')
dev.off()
