
# ### identified treatments and how their profiles change under transformations
# using here lambda=10.5 for comparison with untransformed
# for untransformed data, the T2 was already done in LWafula_HotellingT2_calculation script
# for lambda=10.5 transformation, it's done in the lower version of this script
# lower end there's a working version that can be run as all needed is already saved
# do thz for matched treatments
# ### identified treatments and how their profiles change under transformations
# using here lambda=10.5 for comparison with untransformed

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
matched<-'E:/uhasselt/22yrsem/MThesis/report/results-matricesH2_pvalues etc/matchedactivetreatments'

datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/data/glog transformed datasets by lambda - Copy'
prov.manenoz<-"E:/uhasselt/22yrsem/MThesis/report/dataView/manova plots/proveOfTransformationWork"

# saved untransformed data for treatments corresponding to those in lambda= 10.5
untrmatched_105<-'results_untransformed_lambda10.5.csv'
untrmatched_105<- file.path(matched,untrmatched_105)

untrmatched_105<-read.csv(untrmatched_105,sep=',',header = T)
untrmatched_105<-untrmatched_105[,-1]

# saved transformed data
trmatched_105<-'results_transformed_lambda10.5.csv'
trmatched_105<- file.path(matched,trmatched_105)
trmatched_105<-read.csv(trmatched_105,sep=',',header = T)
trmatched_105<-trmatched_105[,-1]
#

# merge the two results
merged<-merge(untrmatched_105,trmatched_105,by=c('Var1','Var2'))

# minimal values for untransformed and higher values of T2 for the transformed: any difference beyond 250000 considered key

prove.data<-merged[merged[['value.y']]-merged[['value.x']]>=250000,]
prove.data$diff<-prove.data[['value.y']]-prove.data[['value.x']]

# untransformed data
load('data/p1601xxPPS_means_preprocessed.Rdata')

#transformed data

load('E:/uhasselt/22yrsem/MThesis/data/glog transformed datasets by lambda - Copy/p1601xxPPS_means_preprocessedfeat_act_10.5.Rdata')
lambdavals<-10.5
# untransformed data: profiles vs the profiles for features after lambda=10.5 transformation
# names of treatments involved for from low T2 in untransformed to high T2 for transformed
treat.prove<-c()

# random pick
for(i in 1:dim(prove.data)[1]){
  treat.prove[1]<-label_value(prove.data[i,1])
  treat.prove[2]<-label_value(prove.data[i,2])
  
  # now see thz from the untransformed and transformed data
  # limit to only the selected profiles and treatments involved
  datp.prove<-
    datp[ datp[['Active_FractionOfRepl']]>=0.5 & datp[['Treatment']] %in% treat.prove, c("Treatment", attr(datp, "features_mRMR_optN"))]
  
  #replicates 
  datp.prove <- melt(datp.prove, id="Treatment")  # convert to long format
  
  # total counts and indexed counts per cpd [9 reps per cmpd for every feature/max of 45 replicates]
  datp.prove$id=1; datp.prove<-data.table(datp.prove)
  
  datp.prove[, totalreplicates := .N, by = c('Treatment','variable')]
  datp.prove[, indexreplicates := cumsum(id), by = c('Treatment','variable')]
  datp.prove=datp.prove[,-c('id','totalreplicates')]
  
  #renaming
  datp.prove=
    plyr::rename(datp.prove,c('Treatment'='Treatment','variable'='feature','value'='value'))
  
  # rename TRT levels
  datp.prove$Treatment=  gsub("in Colon-CytoSkeleton-Endosomes","",datp.prove$Treatment)
  
  #remove the .zGSLC ending in feature names
  datp.prove$feature=  gsub(".mean.zGSLC","",datp.prove$feature)
  datp.prove$Treatment=as.factor(datp.prove$Treatment)
  
  # untransformed? ; lambda10.5
  glogdatp.prove<-
    glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 & glogdatp[['Treatment']] %in% treat.prove, c("Treatment", attr(glogdatp, "features_mRMR_optN"))]
  
  #replicates per compound
  glogdatp.prove <- melt(glogdatp.prove, id="Treatment")  # convert to long format
  
  # total counts and indexed counts per cpd [9 reps per cmpd for every feature/max of 45 replicates]
  glogdatp.prove$id=1; glogdatp.prove<-data.table(glogdatp.prove)
  
  glogdatp.prove[, totalreplicates := .N, by = c('Treatment','variable')]
  glogdatp.prove[, indexreplicates := cumsum(id), by = c('Treatment','variable')]
  glogdatp.prove=glogdatp.prove[,-c('id','totalreplicates')]
  
  #renaming
  glogdatp.prove=
    plyr::rename(glogdatp.prove,c('Treatment'='Treatment','variable'='feature','value'='value'))
  
  # rename TRT levels
  glogdatp.prove$Treatment=  gsub("in Colon-CytoSkeleton-Endosomes","",glogdatp.prove$Treatment)
  
  #remove the .zGSLC ending in feature names
  glogdatp.prove$feature=  gsub(".zGSLC","",glogdatp.prove$feature)
  glogdatp.prove$Treatment=as.factor(glogdatp.prove$Treatment)
  
  # untransformed plot
  setDT(datp.prove)[, id := .GRP, by = feature]
  puntr<-paste0('graphs/diff250000/lambda',lambdavals,'untr_diff_250000_trt_combn',i,'.pdf')
  puntr<-file.path(prov.manenoz,puntr)
  pdf(puntr)
  par(xpd = T, mar = par()$mar + c(3,0,0,1))
  plot(datp.prove[datp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                  datp.prove$indexreplicates==1,]$id,
       datp.prove[datp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                  datp.prove$indexreplicates==1,]$value,type='n',ylab='replicate value',xlab='',
       xlim=c(1, length(table(datp.prove$feature))),ylim=c(min(datp.prove$value)-1,max(datp.prove$value))+1,xaxt = "n",main='untransformed')
  
  for(k in 1:length(table(datp.prove$indexreplicates))){
    lines(datp.prove[datp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                      datp.prove$indexreplicates==k,]$id,
          datp.prove[datp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                     datp.prove$indexreplicates==k,]$value, col='black')
  }
  for(k in 1:length(table(datp.prove$indexreplicates))){
    lines(datp.prove[datp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[2]) & 
                     datp.prove$indexreplicates==k,]$id,
          datp.prove[datp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[2]) & 
                    datp.prove$indexreplicates==k,]$value, col='red')
  }
  
  axis(1, at=1:length(table(datp.prove$feature)), labels=FALSE)
  text(x=1:length(table(datp.prove$feature)),par("usr")[3]-1,adj = 1,xpd=NA, labels=names(table(datp.prove$feature)), srt=90,
       cex = 0.5)
  legend(3,(max(datp.prove$value)+3),
                    legend=c(gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]),gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[2])),
         col=c('black','red'),bty='n',lty=1:1,ncol = 2)
  #restore mar to default value 
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()

  # transformed plot
  setDT(glogdatp.prove)[, id := .GRP, by = feature]
  ptr<-paste0('graphs/diff250000/lambda',lambdavals,'tr_diff_250000_trt_combn',i,'.pdf')
  ptr<-file.path(prov.manenoz,ptr)
  pdf(ptr)
  par(xpd = T, mar = par()$mar + c(3,0,0,1))
  plot(glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                    glogdatp.prove$indexreplicates==1,]$id,
       glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                    glogdatp.prove$indexreplicates==1,]$value,type='n',ylab='replicate value',xlab='',
       xlim=c(1, length(table(glogdatp.prove$feature))),ylim=c(min(glogdatp.prove$value)-1,max(glogdatp.prove$value))+1,xaxt = "n",
       main=eval(bquote(expression(lambda~'='  ~.(lambdavals)~tranformation))))
  
  for(k in 1:length(table(glogdatp.prove$indexreplicates))){
    lines(glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                       glogdatp.prove$indexreplicates==k,]$id,
          glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                       glogdatp.prove$indexreplicates==k,]$value, col='black')
  }
  for(k in 1:length(table(glogdatp.prove$indexreplicates))){
    lines(glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[2]) & 
                       glogdatp.prove$indexreplicates==k,]$id,
          glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[2]) & 
                       glogdatp.prove$indexreplicates==k,]$value, col='red')
  }
  
  axis(1, at=1:length(table(glogdatp.prove$feature)), labels=FALSE)
  text(x=1:length(table(glogdatp.prove$feature)),par("usr")[3]-1,adj = 1,xpd=NA, labels=names(table(glogdatp.prove$feature)), srt=90,
       cex = 0.5)
  legend(3,(max(glogdatp.prove$value)+3),
         legend=c(gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]),gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[2])),
         col=c('black','red'),bty='n',lty=1:1,ncol = 2)
  #restore mar to default value 
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
  
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #

# those that showed little improvement judging from their diff in T2 pre- and post-tranformation
# ### identified treatments and how their profiles change under transformations
# using here lambda=10.5 for comparison with untransformed
# for untransformed data, the T2 was already done in LWafula_HotellingT2_calculation script
# for lambda=10.5 transformation, it's done in the lower version of this script
# lower end there's a working version that can be run as all needed is already saved
# do thz for matched treatments
# ### identified treatments and how their profiles change under transformations
# using here lambda=10.5 for comparison with untransformed

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
matched<-'E:/uhasselt/22yrsem/MThesis/report/results-matricesH2_pvalues etc/matchedactivetreatments'

datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/data/glog transformed datasets by lambda - Copy'
prov.manenoz<-"E:/uhasselt/22yrsem/MThesis/report/dataView/manova plots/proveOfTransformationWork"

# saved untransformed data for treatments corresponding to those in lambda= 10.5
untrmatched_105<-'results_untransformed_lambda10.5.csv'
untrmatched_105<- file.path(matched,untrmatched_105)

untrmatched_105<-read.csv(untrmatched_105,sep=',',header = T)
untrmatched_105<-untrmatched_105[,-1]

# saved transformed data
trmatched_105<-'results_transformed_lambda10.5.csv'
trmatched_105<- file.path(matched,trmatched_105)
trmatched_105<-read.csv(trmatched_105,sep=',',header = T)
trmatched_105<-trmatched_105[,-1]
#

# merge the two results
merged<-merge(untrmatched_105,trmatched_105,by=c('Var1','Var2'))

# pick out those treatment pairs that showed the least difference b4 and after transformation: keep two that were least different 

prove.data<-merged[abs(merged[['value.y']]-merged[['value.x']])<=0.5,]
prove.data$diff<-prove.data[['value.y']]-prove.data[['value.x']]; prove.data<-data.table(prove.data)
prove.data<-prove.data[order(abs(diff)),]; prove.data<-prove.data[1:4,]
# untransformed data
load('data/p1601xxPPS_means_preprocessed.Rdata')

#transformed data

load('E:/uhasselt/22yrsem/MThesis/data/glog transformed datasets by lambda - Copy/p1601xxPPS_means_preprocessedfeat_act_10.5.Rdata')
lambdavals<-10.5
# untransformed data: profiles vs the profiles for features after lambda=10.5 transformation
# names of treatments involved for from low T2 in untransformed to high T2 for transformed
treat.prove<-c()

# random pick
for(i in 1:dim(prove.data)[1]){
  treat.prove[1]<-label_value(prove.data[i,1])
  treat.prove[2]<-label_value(prove.data[i,2])
  
  # now see thz from the untransformed and transformed data
  # limit to only the selected profiles and treatments involved
  datp.prove<-
    datp[ datp[['Active_FractionOfRepl']]>=0.5 & datp[['Treatment']] %in% treat.prove, c("Treatment", attr(datp, "features_mRMR_optN"))]
  
  #replicates 
  datp.prove <- melt(datp.prove, id="Treatment")  # convert to long format
  
  # total counts and indexed counts per cpd [9 reps per cmpd for every feature/max of 45 replicates]
  datp.prove$id=1; datp.prove<-data.table(datp.prove)
  
  datp.prove[, totalreplicates := .N, by = c('Treatment','variable')]
  datp.prove[, indexreplicates := cumsum(id), by = c('Treatment','variable')]
  datp.prove=datp.prove[,-c('id','totalreplicates')]
  
  #renaming
  datp.prove=
    plyr::rename(datp.prove,c('Treatment'='Treatment','variable'='feature','value'='value'))
  
  # rename TRT levels
  datp.prove$Treatment=  gsub("in Colon-CytoSkeleton-Endosomes","",datp.prove$Treatment)
  
  #remove the .zGSLC ending in feature names
  datp.prove$feature=  gsub(".mean.zGSLC","",datp.prove$feature)
  datp.prove$Treatment=as.factor(datp.prove$Treatment)
  
  # untransformed? ; lambda10.5
  glogdatp.prove<-
    glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 & glogdatp[['Treatment']] %in% treat.prove, c("Treatment", attr(glogdatp, "features_mRMR_optN"))]
  
  #replicates per compound
  glogdatp.prove <- melt(glogdatp.prove, id="Treatment")  # convert to long format
  
  # total counts and indexed counts per cpd [9 reps per cmpd for every feature/max of 45 replicates]
  glogdatp.prove$id=1; glogdatp.prove<-data.table(glogdatp.prove)
  
  glogdatp.prove[, totalreplicates := .N, by = c('Treatment','variable')]
  glogdatp.prove[, indexreplicates := cumsum(id), by = c('Treatment','variable')]
  glogdatp.prove=glogdatp.prove[,-c('id','totalreplicates')]
  
  #renaming
  glogdatp.prove=
    plyr::rename(glogdatp.prove,c('Treatment'='Treatment','variable'='feature','value'='value'))
  
  # rename TRT levels
  glogdatp.prove$Treatment=  gsub("in Colon-CytoSkeleton-Endosomes","",glogdatp.prove$Treatment)
  
  #remove the .zGSLC ending in feature names
  glogdatp.prove$feature=  gsub(".zGSLC","",glogdatp.prove$feature)
  glogdatp.prove$Treatment=as.factor(glogdatp.prove$Treatment)
  
  # untransformed plot
  setDT(datp.prove)[, id := .GRP, by = feature]
  puntr<-paste0('graphs/diffless0.5/lambda',lambdavals,'untr_diff_less0point5_trt_combn',i,'.pdf')
  puntr<-file.path(prov.manenoz,puntr)
  pdf(puntr)
  par(xpd = T, mar = par()$mar + c(3,0,0,1))
  plot(datp.prove[datp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                    datp.prove$indexreplicates==1,]$id,
       datp.prove[datp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                    datp.prove$indexreplicates==1,]$value,type='n',ylab='replicate value',xlab='',
       xlim=c(1, length(table(datp.prove$feature))),ylim=c(min(datp.prove$value)-1,1.25*max(datp.prove$value)),xaxt = "n",main='untransformed')
  
  for(k in 1:length(table(datp.prove$indexreplicates))){
    lines(datp.prove[datp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                       datp.prove$indexreplicates==k,]$id,
          datp.prove[datp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                       datp.prove$indexreplicates==k,]$value, col='black')
  }
  for(k in 1:length(table(datp.prove$indexreplicates))){
    lines(datp.prove[datp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[2]) & 
                       datp.prove$indexreplicates==k,]$id,
          datp.prove[datp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[2]) & 
                       datp.prove$indexreplicates==k,]$value, col='red')
  }
  
  axis(1, at=1:length(table(datp.prove$feature)), labels=FALSE)
  text(x=1:length(table(datp.prove$feature)),par("usr")[3]-1,adj = 1,xpd=NA, labels=names(table(datp.prove$feature)), srt=90,
       cex = 0.5)
  legend(3,(1.248*max(datp.prove$value)),
         legend=c(gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]),gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[2])),
         col=c('black','red'),bty='n',lty=1:1,ncol = 2)
  #restore mar to default value 
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
  
  # transformed plot
  setDT(glogdatp.prove)[, id := .GRP, by = feature]
  ptr<-paste0('graphs/diffless0.5/lambda',lambdavals,'tr_diff_less0point5_trt_combn',i,'.pdf')
  ptr<-file.path(prov.manenoz,ptr)
  pdf(ptr)
  par(xpd = T, mar = par()$mar + c(3,0,0,1))
  plot(glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                        glogdatp.prove$indexreplicates==1,]$id,
       glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                        glogdatp.prove$indexreplicates==1,]$value,type='n',ylab='replicate value',xlab='',
       xlim=c(1, length(table(glogdatp.prove$feature))),ylim=c(min(glogdatp.prove$value)-1,1.25*max(glogdatp.prove$value)),xaxt = "n",
       main=eval(bquote(expression(lambda~'='  ~.(lambdavals)~tranformation))))
  
  for(k in 1:length(table(glogdatp.prove$indexreplicates))){
    lines(glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                           glogdatp.prove$indexreplicates==k,]$id,
          glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]) & 
                           glogdatp.prove$indexreplicates==k,]$value, col='black')
  }
  for(k in 1:length(table(glogdatp.prove$indexreplicates))){
    lines(glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[2]) & 
                           glogdatp.prove$indexreplicates==k,]$id,
          glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[2]) & 
                           glogdatp.prove$indexreplicates==k,]$value, col='red')
  }
  
  axis(1, at=1:length(table(glogdatp.prove$feature)), labels=FALSE)
  text(x=1:length(table(glogdatp.prove$feature)),par("usr")[3]-1,adj = 1,xpd=NA, labels=names(table(glogdatp.prove$feature)), srt=90,
       cex = 0.5)
  legend(3,(max(glogdatp.prove$value)*1.248),
         legend=c(gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[1]),gsub("in Colon-CytoSkeleton-Endosomes","",treat.prove[2])),
         col=c('black','red'),bty='n',lty=1:1,ncol = 2)
  #restore mar to default value 
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
  
}


