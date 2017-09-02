
# ### identified treatments and how their profiles change under transformations
# using here lambda=10.5 for comparison with untransformed

# lower end there's a working version that can be run as all needed is already saved
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
datawdmatrixpvals<-'E:/uhasselt/22yrsem/MThesis/report/results-matricesH2_pvalues etc'

datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/data/glog transformed datasets by lambda - Copy'
prov.manenoz<-"E:/uhasselt/22yrsem/MThesis/report/dataView/manova plots/proveOfTransformationWork"
# 707 treatments, top 10 features
# for every we need 18 obs (18 for each treatment)
load('data/p1601xxPPS_means_preprocessed.Rdata')
datpOpt.feature<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features_mRMR_optN")[1:10])]


# distribution of the Hoteling's T2 statistic and its p-value
treatkp<-names(table(datpOpt.feature$Treatment))
results_trtIDs<-array(NA,dim=c(length(names(table(datpOpt.feature$Treatment))),length(names(table(datpOpt.feature$Treatment)))))


for (j in 1:length(treatkp)) {
  for (i in 1:length(treatkp)) {
    if(j<i) { #you may as well do for (j!=i) though it is just a mirror-image
      fit_ij<-datpOpt.feature[datpOpt.feature$Treatment %in% treatkp[c(j,i)],]
      Y<-as.matrix(cbind(fit_ij[,2:11]),ncol=10,nrow=dim(datpOpt.feature)[1],byrow=T)
      treat_i<-as.factor(fit_ij$Treatment)
      #running MANOVA
      (fit_ij <- manova(Y ~ treat_i))
      colnames(results_trtIDs)<-paste(treatkp,sep = '')
      rownames(results_trtIDs)<-paste(treatkp,sep = '')
      results_trtIDs[j,i]<-summary(fit_ij, test = "Hotelling-Lawley")$stats[1,2]
      print(summary(fit_ij, test = "Hotelling-Lawley")$stats[1,2])
    }
    else {print(paste('i=',i,'j=',j,'already done/no need to do'))}
  }
}

ttrt_namedID<-results_trtIDs
#save


untr<-'data/untransformed.csv'
untr<- file.path(prov.manenoz,untr)
results_trtIDs<-melt(results_trtIDs)
results_trtIDs<-results_trtIDs[complete.cases(results_trtIDs[['value']]),]
results_trtIDs$gr<-'0'

write.csv(results_trtIDs,untr)

#
# into one line for plotting
ttrt_namedID<-melt(ttrt_namedID)
ttrt_namedID<-ttrt_namedID[complete.cases(ttrt_namedID[['value']]),]
ttrt_namedID$gr<-'0'

# lambda 10.5 Hoteling's 
file.names <- dir(datawdtranslambda, pattern ="p1601xxPPS_means_preprocessedfeat_act_")
i=4 #check that this refers to lambda=10.5
  fileb1<-file.path(datawdtranslambda,file.names[i])
  load(fileb1)
  #keep the active treatments, and the top 10 features selected
  datpopt.featurelambda<-glogdatp[glogdatp[['Active_FractionOfRepl']]>0.5 ,c("Treatment", attr(glogdatp, "features_mRMR_optN")[1:10])]
  rm(list=c('glogdatp','fileb1'))
  
  # distribution of the Hoteling's T2 statistic and its p-value
  treatkp<-names(table(datpopt.featurelambda$Treatment)) #703 from the 707 fro untransformed
  #lambda value
  
  lambdaval=   gsub("p1601xxPPS_means_preprocessedfeat_act_","",file.names[i]); lambdaval<-gsub(".RData","",lambdaval); lambdaval
  #resultsname<-paste0('results',lambdaval)
  resultsname<-array(NA,dim=c(length(names(table(datpopt.featurelambda$Treatment))),length(names(table(datpopt.featurelambda$Treatment)))))
  
  for (j in 1:length(treatkp)) {
    for (i in 1:length(treatkp)) {
      if(j<i) { #you may as well do for (j!=i) though it is just a mirror-image
        fit_ij<-datpopt.featurelambda[datpopt.featurelambda$Treatment %in% treatkp[c(j,i)],]
        Y<-as.matrix(cbind(fit_ij[,2:11]),ncol=10,nrow=dim(datpopt.featurelambda)[1],byrow=T)
        treat_i<-as.factor(fit_ij$Treatment)
        #running MANOVA
        (fit_ij <- manova(Y ~ treat_i))
        colnames(resultsname)<-paste(treatkp,sep = '')
        rownames(resultsname)<-paste(treatkp,sep = '') 
        resultsname[j,i]<-summary(fit_ij, test = "Hotelling-Lawley")$stats[1,2]
        
      }
      else {print(paste('i=',i,'j=',j,'already done/no need to do'))}
    }
  }
  
  resultsname10_5<-resultsname
  
  #save
  tr<-'data/transformed_lambda10_5.csv'
  tr<- file.path(prov.manenoz,tr)
  resultsname<-melt(resultsname)
  resultsname<-resultsname[complete.cases(resultsname[['value']]),]
  resultsname$gr<-'10.5'
  write.csv(resultsname,tr)
  
  #
  # into one line for plotting
  resultsname10_5<-melt(resultsname10_5)
  resultsname10_5<-resultsname10_5[complete.cases(resultsname10_5[['value']]),]
  resultsname10_5$gr<-'10.5'
  
# merge the two results
  merged<-merge(ttrt_namedID,resultsname10_5,by=c('Var1','Var2'))
  hi<-paste0('report/dataView/manova plots/individual lambdas/hotellings_untransformed_vs',lambdaval,'.pdf')
  pdf(hi)
plot(merged[['value.x']],merged[['value.y']],pch=19,col='black',cex=0.65,
     xlab=expression(paste('untransformed', ' ', T^2,' ','value')),ylab='')

title(ylab=expression(paste('transformed',' ','(',lambda,' ','= 10.5',')', ' ', T^2,' ','value')),family='',line = 2.5)
dev.off()

# identify the one that led to high changes in T^2 values
# View( merged[merged[['value.y']]==max(merged[['value.y']])|merged[['value.x']]==max(merged[['value.x']]),])

prove.data<-merged[merged[['value.y']]==max(merged[['value.y']])|merged[['value.x']]==max(merged[['value.x']]),]
# untransformed data: profiles vs the profiles for features after lambda=10.5 transformation
# names of treatments involved for from low T2 in untransformed to high T2 for transformed
treat.prove<-c()
# 1stly : from low T2-> max T2
treat.prove[1]<-label_value(prove.data[2,1])
treat.prove[2]<-label_value(prove.data[2,2])

# then from maximal T2 -> minimal T2
treat.prove[3]<-label_value(prove.data[1,1])
treat.prove[4]<-label_value(prove.data[1,2])

# now see thz from the untransformed and transformed data
# limit to only the selected profiles and treatments involved
datp.prove<-
datp[ datp[['Active_FractionOfRepl']]>0.5 & datp[['Treatment']] %in% treat.prove, c("Treatment", attr(datp, "features_mRMR_optN"))]

# profiles b4 and after: 2 of the active compounds
trts=names(table(datp.prove[['Treatment']]))

datp_provemin2max=
  data.table(subset(datp.prove, datp.prove$Treatment %in% trts[3:4]))

#replicates per compound now
datp_provemin2max <- melt(datp_provemin2max, id="Treatment")  # convert to long format

# total counts and indexed counts per cpd [9 reps per cmpd for every feature/max of 45 replicates]
datp_provemin2max$id=1

datp_provemin2max[, totalreplicates := .N, by = c('Treatment','variable')]
datp_provemin2max[, indexreplicates := cumsum(id), by = c('Treatment','variable')]
datp_provemin2max=datp_provemin2max[,-c('id','totalreplicates')]

#renaming
datp_provemin2max=
  plyr::rename(datp_provemin2max,c('Treatment'='Treatment','variable'='feature','value'='value'))

# rename TRT levels
datp_provemin2max$Treatment=  gsub("in Colon-CytoSkeleton-Endosomes","",datp_provemin2max$Treatment)

#remove the .zGSLC ending in feature names
datp_provemin2max$feature=  gsub(".zGSLC","",datp_provemin2max$feature)

#plots of replicates
datp_provemin2max$Treatment=as.factor(datp_provemin2max$Treatment)
hi<-paste0('untrans_minT2_to_lambda',lambdaval,'maxT2.pdf')
hi<-file.path(prov.manenoz,hi)
trans_minT2_to_maxT2=
  ggplot(data = datp_provemin2max,aes(x = feature, y = value,
       group=interaction(indexreplicates,Treatment),col=Treatment)) + geom_point(pch=19,
 cex=0.55) + theme_bw() + geom_line(lwd=0.6)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab('feature')+ylab('replicate value')+ guides(fill=guide_legend(ncol=2)) + 
  theme(legend.position="top",legend.text=element_text(size=10)) 
# ggtitle('features replicate value for 5 compounds in Colon-CytoSkeleton-Endosomes')

ggsave(hi,trans_minT2_to_maxT2)

# after transformation with lambda10.5
load('E:/uhasselt/22yrsem/MThesis/data/glog transformed datasets by lambda - Copy/p1601xxPPS_means_preprocessedfeat_act_10.5.Rdata')
glogdatp.prove<-
  glogdatp[glogdatp[['Active_FractionOfRepl']]>0.5 & glogdatp[['Treatment']] %in% treat.prove, c("Treatment", attr(glogdatp, "features_mRMR_optN"))]

glogdatp_provemin2max=
  data.table(subset(glogdatp.prove, glogdatp.prove$Treatment %in% trts[3:4]))

#replicates per compound now
glogdatp_provemin2max <- melt(glogdatp_provemin2max, id="Treatment")  # convert to long format

# total counts and indexed counts per cpd [9 reps per cmpd for every feature/max of 45 replicates]
glogdatp_provemin2max$id=1

glogdatp_provemin2max[, totalreplicates := .N, by = c('Treatment','variable')]
glogdatp_provemin2max[, indexreplicates := cumsum(id), by = c('Treatment','variable')]
glogdatp_provemin2max=glogdatp_provemin2max[,-c('id','totalreplicates')]

#renaming
glogdatp_provemin2max=
  plyr::rename(glogdatp_provemin2max,c('Treatment'='Treatment','variable'='feature','value'='value'))

# rename TRT levels
glogdatp_provemin2max$Treatment=  gsub("in Colon-CytoSkeleton-Endosomes","",glogdatp_provemin2max$Treatment)

#remove the .zGSLC ending in feature names
glogdatp_provemin2max$feature=  gsub(".zGSLC","",glogdatp_provemin2max$feature)

#plots of replicates
glogdatp_provemin2max$Treatment=as.factor(glogdatp_provemin2max$Treatment)
hi<-paste0('untrans_minT2_to_lambda',lambdaval,'maxT2_transformedversion.pdf')
hi<-file.path(prov.manenoz,hi)
trans_minT2_to_maxT2_transformed=
  ggplot(data = glogdatp_provemin2max,aes(x = feature, y = value,
          group=interaction(indexreplicates,Treatment),col=Treatment)) + geom_point(pch=19,
  cex=0.55) + theme_bw() + geom_line(lwd=0.6)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab('feature')+ylab('replicate value')+ guides(fill=guide_legend(ncol=2)) + 
  theme(legend.position="top",legend.text=element_text(size=10)) 
# ggtitle('features replicate value for 5 compounds in Colon-CytoSkeleton-Endosomes')

ggsave(hi,trans_minT2_to_maxT2_transformed)

# combine them
hicomb<-paste0('untrans_minT2_to_lambda',lambdaval,'maxT2_combined.pdf')
hicomb<-file.path(prov.manenoz,hicomb)

datp_provemin2max$gr<-'before transformation'
glogdatp_provemin2max$gr<-'after transformation'
datacomb<-rbind(datp_provemin2max,glogdatp_provemin2max)

#plot
hicomb0<-
  ggplot(data = datacomb,aes(x = feature, y = value,
group=interaction(indexreplicates,Treatment),col=Treatment)) + geom_point(pch=19,
cex=0.55) + theme_bw() + geom_line(lwd=0.6)+ theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5,color="black"))+
  xlab('feature')+ylab('replicate value')+ guides(fill=guide_legend(ncol=2)) + 
  theme(legend.position="top",legend.text=element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5))+facet_grid(.~gr)+
theme(axis.title=element_text(size=15, color="black",face="bold"))+
theme(panel.background = element_blank())

ggsave(hicomb,hicomb0)

## ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' ###########
# one that changed from maximum under untransformed to minimum

datp_provemax2min=
  data.table(subset(datp.prove, datp.prove$Treatment %in% trts[1:2]))

#replicates per compound now
datp_provemax2min <- melt(datp_provemax2min, id="Treatment")  # convert to long format

# total counts and indexed counts per cpd [9 reps per cmpd for every feature/max of 45 replicates]
datp_provemax2min$id=1

datp_provemax2min[, totalreplicates := .N, by = c('Treatment','variable')]
datp_provemax2min[, indexreplicates := cumsum(id), by = c('Treatment','variable')]
datp_provemax2min=datp_provemax2min[,-c('id','totalreplicates')]

#renaming
datp_provemax2min=
  plyr::rename(datp_provemax2min,c('Treatment'='Treatment','variable'='feature','value'='value'))

# rename TRT levels
datp_provemax2min$Treatment=  gsub("in Colon-CytoSkeleton-Endosomes","",datp_provemax2min$Treatment)

#remove the .zGSLC ending in feature names
datp_provemax2min$feature=  gsub(".zGSLC","",datp_provemax2min$feature)

#plots of replicates
datp_provemax2min$Treatment=as.factor(datp_provemax2min$Treatment)
hi<-paste0('untrans_maxT2_to_lambda',lambdaval,'minT2.pdf')
hi<-file.path(prov.manenoz,hi)
trans_maxT2_to_minT2=
  ggplot(data = datp_provemax2min,aes(x = feature, y = value,
                                      group=interaction(indexreplicates,Treatment),col=Treatment)) + geom_point(pch=19,
                                                                                                                cex=0.55) + theme_bw() + geom_line(lwd=0.6)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab('feature')+ylab('replicate value')+ guides(fill=guide_legend(ncol=2)) + 
  theme(legend.position="top",legend.text=element_text(size=10)) 
# ggtitle('features replicate value for 5 compounds in Colon-CytoSkeleton-Endosomes')

ggsave(hi,trans_maxT2_to_minT2)

# after transformation with lambda10.
glogdatp.prove0<-
  glogdatp[glogdatp[['Active_FractionOfRepl']]>0.5 & glogdatp[['Treatment']] %in% treat.prove, c("Treatment", attr(glogdatp, "features_mRMR_optN"))]

glogdatp_provemax2min=
  data.table(subset(glogdatp.prove0, glogdatp.prove0$Treatment %in% trts[1:2]))

#replicates per compound now
glogdatp_provemax2min <- melt(glogdatp_provemax2min, id="Treatment")  # convert to long format

# total counts and indexed counts per cpd [9 reps per cmpd for every feature/max of 45 replicates]
glogdatp_provemax2min$id=1

glogdatp_provemax2min[, totalreplicates := .N, by = c('Treatment','variable')]
glogdatp_provemax2min[, indexreplicates := cumsum(id), by = c('Treatment','variable')]
glogdatp_provemax2min=glogdatp_provemax2min[,-c('id','totalreplicates')]

#renaming
glogdatp_provemax2min=
  plyr::rename(glogdatp_provemax2min,c('Treatment'='Treatment','variable'='feature','value'='value'))

# rename TRT levels
glogdatp_provemax2min$Treatment=  gsub("in Colon-CytoSkeleton-Endosomes","",glogdatp_provemax2min$Treatment)

#remove the .zGSLC ending in feature names
glogdatp_provemax2min$feature=  gsub(".zGSLC","",glogdatp_provemax2min$feature)

#plots of replicates
glogdatp_provemax2min$Treatment=as.factor(glogdatp_provemax2min$Treatment)
hi<-paste0('untrans_minT2_to_lambda',lambdaval,'maxT2_transformedversion.pdf')
hi<-file.path(prov.manenoz,hi)
trans_minT2_to_maxT2_transformed=
  ggplot(data = glogdatp_provemax2min,aes(x = feature, y = value,
                                          group=interaction(indexreplicates,Treatment),col=Treatment)) + geom_point(pch=19,
                                                                                                                    cex=0.55) + theme_bw() + geom_line(lwd=0.6)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab('feature')+ylab('replicate value')+ guides(fill=guide_legend(ncol=2)) + 
  theme(legend.position="top",legend.text=element_text(size=10)) 
# ggtitle('features replicate value for 5 compounds in Colon-CytoSkeleton-Endosomes')

ggsave(hi,trans_minT2_to_maxT2_transformed)

# combine them
hicomb1.1<-paste0('untrans_maxT2_to_lambda',lambdaval,'minT2_combined.pdf')
hicomb1.1<-file.path(prov.manenoz,hicomb1.1)

datp_provemax2min$gr<-'before transformation'
glogdatp_provemax2min$gr<-'after transformation'
datacomb0<-rbind(datp_provemax2min,glogdatp_provemax2min)

#plot
hicomb1<-
  ggplot(data = datacomb0,aes(x = feature, y = value,group=interaction(indexreplicates,Treatment),col=Treatment)) + 
  geom_point(pch=19, cex=0.55) + theme_bw() + geom_line(lwd=0.6)+ theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5,color="black"))+
  xlab('feature')+ylab('replicate value')+ guides(fill=guide_legend(ncol=2)) + 
  theme(legend.position="top",legend.text=element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5))+facet_grid(.~gr)+
  theme(axis.title=element_text(size=15, color="black",face="bold"))+
  theme(panel.background = element_blank())

ggsave(hicomb1.1,hicomb1)


########## working side - angalia y all features selected with lambda transformation are not means as for untransformed ######

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
datawdmatrixpvals<-'E:/uhasselt/22yrsem/MThesis/report/results-matricesH2_pvalues etc'

datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/data/glog transformed datasets by lambda - Copy'
prov.manenoz<-"E:/uhasselt/22yrsem/MThesis/report/dataView/manova plots/proveOfTransformationWork"

#saved untransformed data
untr<-'data/untransformed.csv'
untr<- file.path(prov.manenoz,untr)

ttrt_namedID<-read.csv(untr,sep=',',header = T)
ttrt_namedID<-ttrt_namedID[,-1]
# saved transformed data
tr<-'data/transformed_lambda10_5.csv'
tr<- file.path(prov.manenoz,tr)
resultsname10_5<-read.csv(tr,sep=',',header = T)
resultsname10_5<-resultsname10_5[,-1]
#

# merge the two results
merged<-merge(ttrt_namedID,resultsname10_5,by=c('Var1','Var2'))

# minimal values for untransformed and higher values of T2 for the transformed: any difference beyond 5000 considered key

prove.data<-merged[abs(merged[['value.y']]-merged[['value.x']])>10000,]
prove.data$diff<-abs(prove.data[['value.y']]-prove.data[['value.x']])

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
    datp[ datp[['Active_FractionOfRepl']]>0.5 & datp[['Treatment']] %in% treat.prove, c("Treatment", attr(datp, "features_mRMR_optN"))]
  
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
  datp.prove$feature=  gsub(".zGSLC","",datp.prove$feature)
  datp.prove$Treatment=as.factor(datp.prove$Treatment)
  # untransformed? ; lambda10.5
  
  glogdatp.prove<-
    glogdatp[glogdatp[['Active_FractionOfRepl']]>0.5 & glogdatp[['Treatment']] %in% treat.prove, c("Treatment", attr(glogdatp, "features_mRMR_optN"))]
  
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
  
  #plot
  datp.prove$gr<-0
  glogdatp.prove$gr<-1
  datacomb<-rbind(datp.prove,glogdatp.prove)
  datacomb$gr <- factor(datacomb$gr, levels=c('0','1'),labels=c('untransformed','transformed'))
  #plot
  nameplot<-paste0('graphs/lambda',lambdavals,'combined_diff_2000_trt_combn',i,'.pdf')
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



############################################################################################################################################

