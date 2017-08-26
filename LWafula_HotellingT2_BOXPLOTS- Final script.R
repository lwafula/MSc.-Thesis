# DATASET2
#Lwafula
# this just uses matrix results from the LWafula-HotellingT2_calculation-final scripts
# for information on the matrices results, run that script
# BOXPLOTS in order to compare with the results from AUC
# 

# 2. includes results from lambdas [9,9.5, and 10.5 T2 statistics compared to untransformed -with related treatments identified]
rm(list=ls())

setwd('E:/uhasselt/22yrsem/MThesis/DATASET2')
library(hexbin)
library(ggplot2)
library(matlab)
library(plyr)
library(dplyr)
library(data.table)
library(flux)
library(WGCNA)
library(knitr) #for kable
library(gplots)

datawd<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data'
datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy'
datawdmatrixpvals<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/results-matricesH2_pvalues etc'

# saved results- untransformed case
results<-read.table('report/results-matricesH2_pvalues etc/July2017results_H2_allTrts10fts_untransformed.csv',header = T,sep=',')
results<-results[,-1]

# into one line for plotting
# results<-melt(results)-already done
results<-results[complete.cases(results[['value']]),]
results$gr<-'0'; results<-results[,c('value','gr')]
#
# add transform distribution lines for Hotelling's
file.names <- dir(datawdmatrixpvals, pattern ="results")
#keep only Hotelling's results for now
file.names<-file.names[!grepl('^(results_)',file.names, ig = TRUE)]; file.names<-file.names[!grepl('^(July2017)',file.names, ig = TRUE)]

lambv<-c()
for (i in 1:length(file.names)){
  h2temp<-file.path(datawdmatrixpvals,file.names[i])
  lambv[i]<-gsub('results','',file.names[i]); lambv[i]<-gsub('.csv','',lambv[i]); lambv[i]
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  
  # into one line for plotting   [h2temp<-melt(h2temp)]
  h2temp<-h2temp[complete.cases(h2temp[['value']]),]
  
  # add transform grp; keep relevant covariates
  h2temp$gr<-lambv[i];  h2temp<-h2temp[,c('value','gr')]
  
  #merge for plotting
  results<-rbind(results,h2temp)
}
rm(h2temp)
#0rdering the lambdas
results$gr<-as.numeric(results$gr)

# save thz data
resultsBoxplots<-'boxpots_data2_data.csv'
resultsBoxplots<-file.path(datawdmatrixpvals,resultsBoxplots)
write.csv(results,resultsBoxplots)

# keep only those with integer values between 0 & 25 since thz others do not seem to change that much
results<-read.table(resultsBoxplots,sep=',',header = T)
results<-results[results$gr<=25,-1]

#plots -all
hi<-paste0('report/dataView/manova plots/individual lambdas/D2_hotellings_allactiveTRTS_Boxplot_all_lambdas','.pdf')
pdf(hi)
plot(value~as.factor(gr),data=results,xlab=expression(paste(lambda,' ','value')),ylab=expression(paste(T^2,' ','value')),pch=19,
     cex=0.25,ylim=c(-7000,max(results$value)))

abline(h=median(results[results$gr==0,]$value),col='red',lwd=2,lty=2)
legend('top',
       legend=c(expression(paste(lambda,' ', '= 0 (untransformed case)')),expression(paste('median for',' ',lambda,' ', '= 0'))),
       pch=c(19,NA),lty=c(NA,2),col=c(NA,'red'),lwd=c(NA,2),bty='n')
dev.off()

# logs
#all
hi<-paste0('report/dataView/manova plots/individual lambdas/D2_loghotellings_allactiveTRTS_Boxplot_all_lambdas','.pdf')
pdf(hi)
plot(log(value)~as.factor(gr),data=results,xlab=expression(paste(lambda,' ','value')),ylab='',pch=19,
     cex=0.25,ylim=c(-log(10),log(max(results$value))), cex.main=0.75,cex.lab=0.85)

abline(h=log(median(results[results$gr==0,]$value)),col='red',lwd=2,lty=2)
legend('bottomleft',
       legend=c(expression(paste(lambda,' ', '= 0 (untransformed case)')),expression(paste('median for',' ',lambda,' ', '= 0'))),
       pch=c(19,NA),lty=c(NA,2),col=c(NA,'red'),lwd=c(NA,2),bty='n')
title(ylab= expression(paste('log','(',T^2,')',' ','value')), line=2.5, cex.lab=1.2, family="")
dev.off()

#limit, anything? [limiting is not that helpful]


