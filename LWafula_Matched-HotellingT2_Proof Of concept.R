
#DATASET2
# ### identified treatments and how their profiles change under transformations
# using here lambda=0.5 for comparison with untransformed
# for untransformed data, the T2 was already done in LWafula_HotellingT2_calculation script
# for lambda=0.5 transformation, it's done in the lower version of this script
# lower end there's a working version that can be run as all needed is already saved
# do thz for matched treatments
# ### identified treatments and how their profiles change under transformations
# using here lambda=0.5 for comparison with untransformed

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

matched<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/results-matricesH2_pvalues etc/matchedactivetreatments'

datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy'
prov.manenoz<-"E:/uhasselt/22yrsem/MThesis/DATASET2/report/dataView/manova plots/proveOfTransformationWork"

# saved untransformed data for treatments corresponding to those in lambda= 0.5
untrmatched_05<-'results_untransformed_lambda0.5.csv'
untrmatched_05<- file.path(matched,untrmatched_05)

untrmatched_05<-read.csv(untrmatched_05,sep=',',header = T)
untrmatched_05<-untrmatched_05[,-1]

# saved transformed data
trmatched_05<-'results_transformed_lambda0.5.csv'
trmatched_05<- file.path(matched,trmatched_05)
trmatched_05<-read.csv(trmatched_05,sep=',',header = T)
trmatched_05<-trmatched_05[,-1]
#

# merge the two results
merged<-merge(untrmatched_05,trmatched_05,by=c('Var1','Var2'))

# minimal values for untransformed and higher values of T2 for the transformed: any difference beyond 200000 considered key

prove.data<-merged[merged[['value.y']]-merged[['value.x']]>=200000,]
prove.data$diff<-prove.data[['value.y']]-prove.data[['value.x']]

# untransformed data
load('E:/uhasselt/22yrsem/MThesis/data/p1602xxPPS_means_preprocessed.Rdata')

#transformed data

load('E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy/p1602xxPPS_means_preprocessedfeat_act_0.5.Rdata')
lambdavals<-0.5
# untransformed data: profiles vs the profiles for features after lambda=0.5 transformation
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
  datp.prove$Treatment=  gsub("in Hepa-ER-Mito","",datp.prove$Treatment)
  
  #remove the .zGSLC ending in feature names
  datp.prove$feature=  gsub(".mean.zGSLC","",datp.prove$feature)
  datp.prove$Treatment=as.factor(datp.prove$Treatment)
  
  # untransformed? ; lambda0.5
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
  glogdatp.prove$Treatment=  gsub("in Hepa-ER-Mito","",glogdatp.prove$Treatment)
  
  #remove the .zGSLC ending in feature names
  glogdatp.prove$feature=  gsub(".zGSLC","",glogdatp.prove$feature)
  glogdatp.prove$Treatment=as.factor(glogdatp.prove$Treatment)
  
  # untransformed plot
  setDT(datp.prove)[, id := .GRP, by = feature]
  puntr<-paste0('graphs/diff200000/D2_lambda',lambdavals,'untr_diff_200000_trt_combn',i,'.pdf')
  puntr<-file.path(prov.manenoz,puntr)
  pdf(puntr)
  par(xpd = T, mar = par()$mar + c(4,0,0,1))
  plot(datp.prove[datp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                    datp.prove$indexreplicates==1,]$id,
       datp.prove[datp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                    datp.prove$indexreplicates==1,]$value,type='n',ylab='replicate value',xlab='',
       xlim=c(1, length(table(datp.prove$feature))),ylim=c(min(datp.prove$value)-1,max(datp.prove$value))+1,xaxt = "n",main='untransformed')
  
  for(k in 1:length(table(datp.prove$indexreplicates))){
    lines(datp.prove[datp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                       datp.prove$indexreplicates==k,]$id,
          datp.prove[datp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                       datp.prove$indexreplicates==k,]$value, col='black')
  }
  for(k in 1:length(table(datp.prove$indexreplicates))){
    lines(datp.prove[datp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[2]) & 
                       datp.prove$indexreplicates==k,]$id,
          datp.prove[datp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[2]) & 
                       datp.prove$indexreplicates==k,]$value, col='red')
  }
  
  axis(1, at=1:length(table(datp.prove$feature)), labels=FALSE)
  text(x=1:length(table(datp.prove$feature)),par("usr")[3]-1,adj = 1,xpd=NA, labels=names(table(datp.prove$feature)), srt=90,
       cex = 0.5)
  legend(3,(max(datp.prove$value)+3),
         legend=c(gsub("in Hepa-ER-Mito","",treat.prove[1]),gsub("in Hepa-ER-Mito","",treat.prove[2])),
         col=c('black','red'),bty='n',lty=1:1,ncol = 2)
  #restore mar to default value 
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
  
  # transformed plot
  setDT(glogdatp.prove)[, id := .GRP, by = feature]
  ptr<-paste0('graphs/diff200000/D2_lambda',lambdavals,'tr_diff_200000_trt_combn',i,'.pdf')
  ptr<-file.path(prov.manenoz,ptr)
  pdf(ptr)
  par(xpd = T, mar = par()$mar + c(4,0,0,1))
  plot(glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                        glogdatp.prove$indexreplicates==1,]$id,
       glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                        glogdatp.prove$indexreplicates==1,]$value,type='n',ylab='replicate value',xlab='',
       xlim=c(1, length(table(glogdatp.prove$feature))),ylim=c(min(glogdatp.prove$value)-1,max(glogdatp.prove$value))+1,xaxt = "n",
       main=eval(bquote(expression(lambda~'='  ~.(lambdavals)~tranformation))))
  
  for(k in 1:length(table(glogdatp.prove$indexreplicates))){
    lines(glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                           glogdatp.prove$indexreplicates==k,]$id,
          glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                           glogdatp.prove$indexreplicates==k,]$value, col='black')
  }
  for(k in 1:length(table(glogdatp.prove$indexreplicates))){
    lines(glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[2]) & 
                           glogdatp.prove$indexreplicates==k,]$id,
          glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[2]) & 
                           glogdatp.prove$indexreplicates==k,]$value, col='red')
  }
  
  axis(1, at=1:length(table(glogdatp.prove$feature)), labels=FALSE)
  text(x=1:length(table(glogdatp.prove$feature)),par("usr")[3]-1,adj = 1,xpd=NA, labels=names(table(glogdatp.prove$feature)), srt=90,
       cex = 0.5)
  legend(3,(max(glogdatp.prove$value)+3),
         legend=c(gsub("in Hepa-ER-Mito","",treat.prove[1]),gsub("in Hepa-ER-Mito","",treat.prove[2])),
         col=c('black','red'),bty='n',lty=1:1,ncol = 2)
  #restore mar to default value 
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
  
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ THOSE WITH VERY SMALL CHANGES IN T2 after and b4 transformation @@@@@@@@@@@@@@@@@@@@@@@ ##

# those that showed little improvement even after tranformation

# ### identified treatments and how their profiles change under transformations
# using here lambda=0.5 for comparison with untransformed
# for untransformed data, the T2 was already done in LWafula_HotellingT2_calculation script
# for lambda=0.5 transformation, it's done in the lower version of this script
# lower end there's a working version that can be run as all needed is already saved
# do thz for matched treatments
# ### identified treatments and how their profiles change under transformations
# using here lambda=0.5 for comparison with untransformed

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

matched<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/results-matricesH2_pvalues etc/matchedactivetreatments'

datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy'
prov.manenoz<-"E:/uhasselt/22yrsem/MThesis/DATASET2/report/dataView/manova plots/proveOfTransformationWork"

# saved untransformed data for treatments corresponding to those in lambda= 0.5
untrmatched_05<-'results_untransformed_lambda0.5.csv'
untrmatched_05<- file.path(matched,untrmatched_05)

untrmatched_05<-read.csv(untrmatched_05,sep=',',header = T)
untrmatched_05<-untrmatched_05[,-1]

# saved transformed data
trmatched_05<-'results_transformed_lambda0.5.csv'
trmatched_05<- file.path(matched,trmatched_05)
trmatched_05<-read.csv(trmatched_05,sep=',',header = T)
trmatched_05<-trmatched_05[,-1]
#

# merge the two results
merged<-merge(untrmatched_05,trmatched_05,by=c('Var1','Var2'))

# those that showed a difference of less than 0.5

prove.data<-data.table(merged[abs(merged[['value.y']]-merged[['value.x']])<=0.5,])
prove.data$diff<-prove.data[['value.y']]-prove.data[['value.x']]; prove.data<-prove.data[order(abs(diff)),] 
prove.data<-prove.data[1:4,]

# untransformed data
load('E:/uhasselt/22yrsem/MThesis/data/p1602xxPPS_means_preprocessed.Rdata')

# # matrix for analysis
# # datp[c(1:4,1530:1534,10000:10004,14684:14688),1:5]
# datp1<-datp[c(1:4,1530:1534,10000:10004,14684:14688),c("Treatment","WELL_ID", attr(datp,"features"))] 
# datp1<-data.table(datp1)
# # remove .mean.ZGSLC
# names(datp1) <- gsub(".mean.zGSLC", " ", names(datp1))
# names(datp1) <- gsub(".zGSLC", " ", names(datp1))
# datp1$Treatment <- gsub(" in Hepa-ER-Mito", " ", datp1$Treatment)
# write.csv(datp1, file="E:/uhasselt/22yrsem/MThesis/report/dataView/EDA/dataMatrix.csv")

#transformed data

load('E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy/p1602xxPPS_means_preprocessedfeat_act_0.5.Rdata')
lambdavals<-0.5
# untransformed data: profiles vs the profiles for features after lambda=0.5 transformation
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
  datp.prove$Treatment=  gsub("in Hepa-ER-Mito","",datp.prove$Treatment)
  
  #remove the .zGSLC ending in feature names
  datp.prove$feature=  gsub(".mean.zGSLC","",datp.prove$feature)
  datp.prove$Treatment=as.factor(datp.prove$Treatment)
  
  # untransformed? ; lambda0.5
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
  glogdatp.prove$Treatment=  gsub("in Hepa-ER-Mito","",glogdatp.prove$Treatment)
  
  #remove the .zGSLC ending in feature names
  glogdatp.prove$feature=  gsub(".zGSLC","",glogdatp.prove$feature)
  glogdatp.prove$Treatment=as.factor(glogdatp.prove$Treatment)
  
  # untransformed plot
  setDT(datp.prove)[, id := .GRP, by = feature]
  puntr<-paste0('graphs/diffless0.5/D2_lambda',lambdavals,'untr_diff_less0point5_trt_combn',i,'.pdf')
  puntr<-file.path(prov.manenoz,puntr)
  pdf(puntr)
  par(xpd = T, mar = par()$mar + c(4,0,0,1))
  plot(datp.prove[datp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                    datp.prove$indexreplicates==1,]$id,
       datp.prove[datp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                    datp.prove$indexreplicates==1,]$value,type='n',ylab='replicate value',xlab='',
       xlim=c(1, length(table(datp.prove$feature))),ylim=c(min(datp.prove$value)-1,max(datp.prove$value))*1.25,
       xaxt = "n",main='untransformed')
  
  for(k in 1:length(table(datp.prove$indexreplicates))){
    lines(datp.prove[datp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                       datp.prove$indexreplicates==k,]$id,
          datp.prove[datp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                       datp.prove$indexreplicates==k,]$value, col='black')
  }
  for(k in 1:length(table(datp.prove$indexreplicates))){
    lines(datp.prove[datp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[2]) & 
                       datp.prove$indexreplicates==k,]$id,
          datp.prove[datp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[2]) & 
                       datp.prove$indexreplicates==k,]$value, col='red')
  }
  
  axis(1, at=1:length(table(datp.prove$feature)), labels=FALSE)
  text(x=1:length(table(datp.prove$feature)),par("usr")[3]-1,adj = 1,xpd=NA, labels=names(table(datp.prove$feature)), srt=90,
       cex = 0.5)
  legend(3,(max(datp.prove$value)*1.248),
         legend=c(gsub("in Hepa-ER-Mito","",treat.prove[1]),gsub("in Hepa-ER-Mito","",treat.prove[2])),
         col=c('black','red'),bty='n',lty=1:1,ncol = 2)
  #restore mar to default value 
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
  
  # transformed plot
  setDT(glogdatp.prove)[, id := .GRP, by = feature]
  ptr<-paste0('graphs/diffless0.5/D2_lambda',lambdavals,'tr_diff_less0point5_trt_combn',i,'.pdf')
  ptr<-file.path(prov.manenoz,ptr)
  pdf(ptr)
  par(xpd = T, mar = par()$mar + c(4,0,0,1))
  plot(glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                        glogdatp.prove$indexreplicates==1,]$id,
       glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                        glogdatp.prove$indexreplicates==1,]$value,type='n',ylab='replicate value',xlab='',
       xlim=c(1, length(table(glogdatp.prove$feature))),ylim=c(min(glogdatp.prove$value)-1,max(glogdatp.prove$value))*1.25,
       xaxt = "n", main=eval(bquote(expression(lambda~'='  ~.(lambdavals)~tranformation))))
  
  for(k in 1:length(table(glogdatp.prove$indexreplicates))){
    lines(glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                           glogdatp.prove$indexreplicates==k,]$id,
          glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[1]) & 
                           glogdatp.prove$indexreplicates==k,]$value, col='black')
  }
  for(k in 1:length(table(glogdatp.prove$indexreplicates))){
    lines(glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[2]) & 
                           glogdatp.prove$indexreplicates==k,]$id,
          glogdatp.prove[glogdatp.prove$Treatment %in% gsub("in Hepa-ER-Mito","",treat.prove[2]) & 
                           glogdatp.prove$indexreplicates==k,]$value, col='red')
  }
  
  axis(1, at=1:length(table(glogdatp.prove$feature)), labels=FALSE)
  text(x=1:length(table(glogdatp.prove$feature)),par("usr")[3]-1,adj = 1,xpd=NA, labels=names(table(glogdatp.prove$feature)), srt=90,
       cex = 0.5)
  legend(3,(max(glogdatp.prove$value)*1.248),
         legend=c(gsub("in Hepa-ER-Mito","",treat.prove[1]),gsub("in Hepa-ER-Mito","",treat.prove[2])),
         col=c('black','red'),bty='n',lty=1:1,ncol = 2)
  #restore mar to default value 
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
  
}



