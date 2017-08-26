
#DATASET2
# ### identified treatments and how their profiles change under transformations
# using here lambda=10.5 for comparison with untransformed
# for untransformed data, the T2 was already done in LWafula_HotellingT2_calculation script
# for lambda=10.5 transformation, it's done in the lower version of this script
# lower end there's a working version that can be run as all needed is already saved


# ### identified treatments and how their profiles change under transformations
# using here lambda=10.5 for comparison with untransformed

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

datawd<-'E:/uhasselt/22yrsem/MThesis/data'
datawdmatrixpvals<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/results-matricesH2_pvalues etc'

datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy'
prov.manenoz<-"E:/uhasselt/22yrsem/MThesis/DATASET2/report/dataView/manova plots/proveOfTransformationWork"

# saved untransformed data
untr<-'data/July2017results_H2_allTrts10fts_untransformed.csv'
untr<- file.path(prov.manenoz,untr)

ttrt_namedID<-read.csv(untr,sep=',',header = T)
ttrt_namedID<-ttrt_namedID[,-1]

# saved transformed data
tr<-'data/results0.5.csv'
tr<- file.path(prov.manenoz,tr)
resultsname0.5<-read.csv(tr,sep=',',header = T)
resultsname0.5<-resultsname0.5[,-1]
#

# merge the two results
merged<-merge(ttrt_namedID,resultsname0.5,by=c('Var1','Var2'))

# minimal values for untransformed and higher values of T2 for the transformed: any difference beyond 250000 considered key

prove.data<-merged[merged[['value.y']]-merged[['value.x']]>=200000,]
prove.data$diff<-prove.data[['value.y']]-prove.data[['value.x']]

# untransformed data
load('E:/uhasselt/22yrsem/MThesis/data/p1602xxPPS_means_preprocessed.Rdata')

#transformed data

load('E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy/p1602xxPPS_means_preprocessedfeat_act_0.5.Rdata')
lambdavals<-0.5
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
  datp.prove$Treatment=  gsub("in Hepa-ER-Mito","",datp.prove$Treatment)
  
  #remove the .zGSLC ending in feature names
  datp.prove$feature=  gsub(".zGSLC","",datp.prove$feature)
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
  glogdatp.prove$Treatment=  gsub("in Hepa-ER-Mito","",glogdatp.prove$Treatment)
  
  #remove the .zGSLC ending in feature names
  glogdatp.prove$feature=  gsub(".zGSLC","",glogdatp.prove$feature)
  glogdatp.prove$Treatment=as.factor(glogdatp.prove$Treatment)
  
  #plot
  datp.prove$gr<-0
  glogdatp.prove$gr<-1
  datacomb<-rbind(datp.prove,glogdatp.prove)
  datacomb$gr <- factor(datacomb$gr, levels=c('0','1'),labels=c('untransformed','transformed'))
  #plot
  nameplot<-paste0('graphs/diff200000/D2_lambda',lambdavals,'combined_diff_200000_trt_combn',i,'.pdf')
  nameplot<-file.path(prov.manenoz,nameplot)
  plots<-ggplot(data = datacomb,aes(x = feature, y = value,group=interaction(indexreplicates,Treatment),col=Treatment)) + 
    geom_point(pch=19,cex=0.55) + theme_bw() + geom_line(lwd=0.6)+ theme(axis.text.x = element_text(angle = 90, hjust = 1,size=3,color="black"))+
    xlab('feature')+ylab('replicate value')+ guides(fill=guide_legend(ncol=2)) + 
    theme(legend.position="top",legend.text=element_text(size=10))+
    theme(plot.title = element_text(hjust = 0.5))+facet_grid(~gr,scales="free", space = 'free')+facet_wrap(~gr, scales = "free_y")+
    theme(axis.title=element_text(size=15, color="black",face="bold"))+
    theme(panel.background = element_blank())
  ggsave(nameplot,plots)
}


# ```````````` for lambda=0.1

#DATASET2
# ### identified treatments and how their profiles change under transformations
# using here lambda=10.5 for comparison with untransformed
# for untransformed data, the T2 was already done in LWafula_HotellingT2_calculation script
# for lambda=10.5 transformation, it's done in the lower version of this script
# lower end there's a working version that can be run as all needed is already saved


# ### identified treatments and how their profiles change under transformations
# using here lambda=10.5 for comparison with untransformed

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

datawd<-'E:/uhasselt/22yrsem/MThesis/data'
datawdmatrixpvals<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/results-matricesH2_pvalues etc'

datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy'
prov.manenoz<-"E:/uhasselt/22yrsem/MThesis/DATASET2/report/dataView/manova plots/proveOfTransformationWork"

# saved untransformed data
untr<-'data/July2017results_H2_allTrts10fts_untransformed.csv'
untr<- file.path(prov.manenoz,untr)

ttrt_namedID<-read.csv(untr,sep=',',header = T)
ttrt_namedID<-ttrt_namedID[,-1]

# saved transformed data
tr<-'data/results0.1.csv'
tr<- file.path(prov.manenoz,tr)
resultsname0.1<-read.csv(tr,sep=',',header = T)
resultsname0.1<-resultsname0.1[,-1]
#

# merge the two results
merged<-merge(ttrt_namedID,resultsname0.1,by=c('Var1','Var2'))

# minimal values for untransformed and higher values of T2 for the transformed: any difference beyond 250000 considered key

prove.data<-merged[merged[['value.y']]-merged[['value.x']]>=200000,]
prove.data$diff<-prove.data[['value.y']]-prove.data[['value.x']]

# untransformed data
load('E:/uhasselt/22yrsem/MThesis/data/p1602xxPPS_means_preprocessed.Rdata')

#transformed data

load('E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy/p1602xxPPS_means_preprocessedfeat_act_0.1.Rdata')
lambdavals<-0.1
# untransformed data: profiles vs the profiles for features after lambda=10.1 transformation
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
  datp.prove$feature=  gsub(".zGSLC","",datp.prove$feature)
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
  glogdatp.prove$Treatment=  gsub("in Hepa-ER-Mito","",glogdatp.prove$Treatment)
  
  #remove the .zGSLC ending in feature names
  glogdatp.prove$feature=  gsub(".zGSLC","",glogdatp.prove$feature)
  glogdatp.prove$Treatment=as.factor(glogdatp.prove$Treatment)
  
  #plot
  datp.prove$gr<-0
  glogdatp.prove$gr<-1
  datacomb<-rbind(datp.prove,glogdatp.prove)
  datacomb$gr <- factor(datacomb$gr, levels=c('0','1'),labels=c('untransformed','transformed'))
  #plot
  nameplot<-paste0('graphs/diff200000/D2_lambda',lambdavals,'combined_diff_200000_trt_combn',i,'.pdf')
  nameplot<-file.path(prov.manenoz,nameplot)
  plots<-ggplot(data = datacomb,aes(x = feature, y = value,group=interaction(indexreplicates,Treatment),col=Treatment)) + 
    geom_point(pch=19,cex=0.55) + theme_bw() + geom_line(lwd=0.6)+ theme(axis.text.x = element_text(angle = 90, hjust = 1,size=3,color="black"))+
    xlab('feature')+ylab('replicate value')+ guides(fill=guide_legend(ncol=2)) + 
    theme(legend.position="top",legend.text=element_text(size=10))+
    theme(plot.title = element_text(hjust = 0.5))+facet_grid(~gr,scales="free", space = 'free')+facet_wrap(~gr, scales = "free_y")+
    theme(axis.title=element_text(size=15, color="black",face="bold"))+
    theme(panel.background = element_blank())
  ggsave(nameplot,plots)
}

