
# aggregated
#Lwafula
# this just uses matrix results from the LWafula-HotellingT2_calculation-final scripts
# for information on the matrices results, run that script
# BOXPLOTS in order to compare with the results from AUC
# 

# 2. includes results from lambdas [9,9.5, and 10.5 T2 statistics compared to untransformed -with related treatments identified]
rm(list=ls())

setwd('E:/uhasselt/22yrsem/MThesis')
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

datawd<-'E:/uhasselt/22yrsem/MThesis/data'
datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/data/glog transformed datasets by lambda - Copy'
datawdmatrixpvals<-'E:/uhasselt/22yrsem/MThesis/report/results-matricesH2_pvalues etc'

# add transform distribution lines for Hotelling's for each of the lambdas [both cell-level and well-level transforms]

lambv<-c()
lambdas=c(0.5,8.5,10,10.5,11,11.5)

for (i in 1:length(lambdas)){
  # saved results- untransformed case
  results<-read.table('report/results-matricesH2_pvalues etc/July2017results_H2_allTrts10fts_untransformed.csv',header = T,sep=',')
  results<-results[,-1]
  
  # into one line for plotting
  # results<-melt(results)-already done
  results<-results[complete.cases(results[['value']]),]
  results$gr<-'0'; results<-results[,c('value','gr')]
  
  #single cell
  h2temp<-file.path(datawdmatrixpvals,paste0('results',lambdas[i],'.csv'))
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  # into one line for plotting
  h2temp<-melt(h2temp)
  h2temp<-h2temp[complete.cases(h2temp[['value']]),]
  
  # add transform grp; keep relevant covariates
  h2temp$gr<-'1';  h2temp<-h2temp[,c('value','gr')]
  
  #merge for plotting
  results<-rbind(results,h2temp)
# well level
  h2temp<-file.path(datawdmatrixpvals,paste0('aggregated_results',lambdas[i],'.csv'))
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  h2temp$gr<-'2';  h2temp<-h2temp[,c('value','gr')]
  
  #merge for plotting
  results<-rbind(results,h2temp)
  rm(h2temp)
  #0rdering the lambdas
  results$gr<-as.numeric(results$gr)
  results$gr<-factor(results$gr,levels=c(0,1,2),labels=c('untransformed','cell level','well level'))
  # logs
  #all
  hi<-paste0('report/dataView/manova plots/individual lambdas/aggregated_loghotellings_allactiveTRTS_Boxplot',lambdas[i],'.pdf')
  pdf(hi)
  plot(log(value)~gr,data=results,ylab='',pch=19,
       cex=0.25,ylim=c(-log(10),log(max(results$value))), cex.main=0.75,cex.lab=0.85)
  
  abline(h=log(median(results[results$gr=='untransformed',]$value)),col='red',lwd=2,lty=2)
  legend('bottomleft',legend=c(expression(paste('median for untransformed'))),lty=c(2),col=c('red'),lwd=c(2),bty='n')
  title(ylab= expression(paste('log','(',T^2,')',' ','value')), line=2.5, cex.lab=1.2, family="")
  dev.off()
  
  }

#limit, anything? [limiting is not that helpful]


