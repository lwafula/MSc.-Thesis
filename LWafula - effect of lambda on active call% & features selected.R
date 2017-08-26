
# LWAFULA

# effect of lambda on the active calls/ number of features selected
# 10-06-2017

# the data used here is from the script named : LWafula- general lambda tranformed  process from startto end- how many minutes.R

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

datawd<-'E:/uhasselt/22yrsem/MThesis/data/glog transformed datasets by lambda'
dataglog<-'E:/uhasselt/22yrsem/MThesis/data/glog_estimates'

#data on the Nfeatures and active call % after alpha transformation
# #save
lambdafeatactv_alphacorrected<-'lambdafeatactv_alphacorrected.csv'
lambdafeatactv_alphacorrected<-file.path(datawd,lambdafeatactv_alphacorrected)
lambdafeatactv_alphacorrected<-read.csv(lambdafeatactv_alphacorrected,header=T,sep=',')
lambdafeatactv_alphacorrected<-lambdafeatactv_alphacorrected[,-1]
lambdafeatactv_alphacorrected<-
  plyr::rename(lambdafeatactv_alphacorrected,c('V1'='lambda','V2'='nfeature','V3'='active_call'))

# #plots
lambdavsoptfeature_alpha<-'report/dataView/lambdavsoptfeature_alpha.pdf'
pdf(lambdavsoptfeature_alpha)
plot(lambdafeatactv_alphacorrected[,1],lambdafeatactv_alphacorrected[,2],pch=19,type='l',xlab=expression(lambda),
     ylab='Optimal features',ylim=c(10,30),lwd=1,lty=1)
abline(h=19,lty=2,col='black',lwd=2)
legend(5.5,30,legend=c(expression(paste('transformed by',' ',lambda,' ','= x-value')),'untransformed'),lty=c(1,2),col=c('black','black'),
       lwd=c(1,2),bty='n')
dev.off()
# 
lambdavsactcall_alpha<-'report/dataView/lambdavsactcall_alpha.pdf'
pdf(lambdavsactcall_alpha)
plot(lambdafeatactv_alphacorrected[,1],lambdafeatactv_alphacorrected[,3],pch=19,type='l',xlab=expression(lambda),
     ylab='active calling %',ylim=c(55,59),lwd=1,lty=1)
abline(h=56.42458,lty=2,col='black',lwd=2)
legend(5.5,59,legend=c(expression(paste('transformed by',' ',lambda,' ','= x-value')),'untransformed'),
       lty=c(1,2),col=c('black','black'),lwd=c(1,2),bty='n')
dev.off()
# 

