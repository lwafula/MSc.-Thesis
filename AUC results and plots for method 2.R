
rm(list=ls())
# UPDATED on 27-07-2017 (LWafula)
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
datp<-
 datp[datp[['Active_FractionOfRepl']]>=0.5,]

#
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
untrans<-'untransformed.pdf' #all UAC values- data.frame(FSres$aucAll)
untrans<-file.path(datawdauc,untrans)
pdf(untrans)
plot(x = FSres)
dev.off()

#distributions for the untransformed cases (replicates vs non-replicates)
#replicates vs non-replicates
for (j in 1:10){
  plotname<-paste0('replicatesvsnonreplicates_run_',j,'_for_untransformed','.pdf')
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
mpya$gr<-'0'
# optN <- FSres$optimalNFeatures # optimal number of features
mpya<-mpya[,c('value','nFeatures','run','gr')]

######### 
#for lambda 10.5 (i=4,); here we do not limit the number of optimal features, just working with active treatments
lambv<-c()
file.names <- dir(datawdtranslambda, pattern ="p1601xxPPS_means_preprocessedfeat_act_")
for (i in 1:length(file.names)){
  fileb1<-file.path(datawdtranslambda,file.names[i])
  lambv[i]<-gsub('p1601xxPPS_means_preprocessedfeat_act_','',file.names[i]); lambv[i]<-gsub('.RData','',lambv[i]); lambv[i]
  
  load(fileb1)
  datpopt.featurelambda<-
    glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,]
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

mpya.1$gr<-lambv[i]
mpya.1<-mpya.1[,c('value','nFeatures','run','gr')]

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
plotname<-paste0('auc_for_lambda',lambv[i],'.pdf')
plotname<-file.path(datawdauc,plotname)
pdf(plotname)
plot(x = FSres.featurelambda, type = "aucVsFeatures")
dev.off()

#replicates vs non-replicates
for (j in 1:10){
  plotname<-paste0('replicatesvsnonreplicates_run_',j,'_for_lambda',lambv[i],'.pdf')
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
write.csv(mpya,file='E:/uhasselt/22yrsem/MThesis/report/dataView/aucmethod/calculated_results/allauc.csv')
#drop for aggregated for now
mpya<-mpya[grep("aggregated_", mpya$gr, invert = TRUE) , ]
lambv<-lambv[grep("aggregated_", lambv, invert = TRUE)]

## combine the plots for all AUC for all lambdas against untransformed
mpya$gr<-as.numeric(mpya$gr)
auc<-'auc.pdf'
auc<-file.path(datawdauc,auc)
pdf(auc)
plot(value~as.factor(gr),data=mpya,xlab=expression(paste(lambda,' ','value')),
ylab='AUC of replicates vs non-replicates across runs',ylim=c(0.97,0.99),main='Evolution of AUC assessing replicability Vs
transformation parameter',pch=19)
abline(h=median(FSres$aucAll),col='red',lwd=2,lty=2)
legend('bottomleft',
       legend=c(expression(paste(lambda,' ', '= 0 (untransformed case)')),expression(paste('median for',' ',lambda,' ', '= 0'))),
       pch=c(19,NA),lty=c(NA,2),col=c(NA,'red'),lwd=c(NA,2),bty='n')
dev.off()

#summaries per group
lambv<-as.numeric(lambv)
for(i in 1:length(lambv)){
summary(mpya[mpya$gr==lambv[i],]$value)
}
