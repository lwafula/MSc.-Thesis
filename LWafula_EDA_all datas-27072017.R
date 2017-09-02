
# `````````````` OBJ.1 : showing the differentiating ability of the selected features ``````````````````````````````````````````
# LWafula 27072017
# Thz script is primarily for all EDA for all data-performed on all the cellLines
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
library(Hotelling) # for Hotelling calculation

datawd<-'E:/uhasselt/22yrsem/MThesis/data'
eda<-'E:/uhasselt/22yrsem/MThesis/report/dataView/EDA'

# how many features are captured at cell-level
plate67070<-'phaedra_singlecell_p160118PPS008__plate67070.Rdata'

plate67070<-file.path(datawd,plate67070)
load(plate67070)
plate67070<-dc

# single cell data 
plate67070 <- data.table(plate67070) 

#Identify the feature columns and remove troublesome features
cn <- names(plate67070)

#grepl returns TRUE if a string contains the pattern, otherwise FALSE [462 features]
features <- cn[grepl('^(Cell_|Cyto_|Nuc_|Ratio_|CellCount|CellConfluency)',
                     cn, ig = TRUE)]

rm(plate67070,dc,features,cn)

# 707 treatments, optimal features [untransformed data]
load('data/p1601xxPPS_means_preprocessed.Rdata')
datpOpt.feature<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features_mRMR_optN"))]
rm(list=c('datp'))

# plot profiles for treatments that are highly similar via Hotelling's T2 statistic+those that are very dissimilar using T2
# idea from the Lwafula_HotellingT2_calculation-Final script
## [cmpd193 @ 3 uM in Colon-CytoSkeleton-Endosomes & cmpd240 @ 9 uM in Colon-CytoSkeleton-Endosomes with T2=322017.5]- large diff
# [cmpd60 @ 1 uM in Colon-CytoSkeleton-Endosomes & cmpd60 @ 3 uM in Colon-CytoSkeleton-Endosomes with T2=7.7]- little diff

trtkp<-c('cmpd193 @ 3 uM in Colon-CytoSkeleton-Endosomes','cmpd240 @ 9 uM in Colon-CytoSkeleton-Endosomes',
         'cmpd60 @ 1 uM in Colon-CytoSkeleton-Endosomes','cmpd60 @ 3 uM in Colon-CytoSkeleton-Endosomes')

# keep only the target ones
datpOpt.feature<-datpOpt.feature[datpOpt.feature[['Treatment']] %in% trtkp,]
table(datpOpt.feature$Treatment)


# maximal difference combination
datpOpt.feature.max<-
  datpOpt.feature[datpOpt.feature[['Treatment']] %in% c('cmpd193 @ 3 uM in Colon-CytoSkeleton-Endosomes',
                                                        'cmpd240 @ 9 uM in Colon-CytoSkeleton-Endosomes'),]

datpOpt.feature.max.long <- melt(datpOpt.feature.max, id="Treatment")  # convert to long format
rm(datpOpt.feature.max)
# total counts and indexed counts per cpd [9 reps per cmpd for every feature]
datpOpt.feature.max.long$id=1; datpOpt.feature.max.long<-data.table(datpOpt.feature.max.long)

datpOpt.feature.max.long[, totalreplicates := .N, by = c('Treatment','variable')]
datpOpt.feature.max.long[, indexreplicates := cumsum(id), by = c('Treatment','variable')]
datpOpt.feature.max.long=datpOpt.feature.max.long[,-c('id','totalreplicates')]

#renaming
datpOpt.feature.max.long=
  plyr::rename(datpOpt.feature.max.long,c('Treatment'='Treatment','variable'='feature','value'='value'))

# rename TRT levels
datpOpt.feature.max.long$Treatment=  gsub("in Colon-CytoSkeleton-Endosomes","",datpOpt.feature.max.long$Treatment)

#plots of replicates
datpOpt.feature.max.long$Treatment=as.factor(datpOpt.feature.max.long$Treatment)
replicatesall_featselected=  ggplot(data = datpOpt.feature.max.long,aes(x = feature, y = value,group=interaction(indexreplicates,Treatment),col=Treatment)) +
  geom_point(pch=19, cex=0.55) + theme_bw() + geom_line(lwd=0.6)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab('feature')+
  ylab('replicate value')+ guides(fill=guide_legend(ncol=2)) + theme(legend.position="top",legend.text=element_text(size=10)) 

ggsave('report/dataView/EDA/EDA_maximal_diffTRT.pdf',replicatesall_featselected)

# minimal diff
datpOpt.feature.min<-
  datpOpt.feature[datpOpt.feature[['Treatment']] %in% c('cmpd60 @ 1 uM in Colon-CytoSkeleton-Endosomes','cmpd60 @ 3 uM in Colon-CytoSkeleton-Endosomes'),]

datpOpt.feature.min.long <- melt(datpOpt.feature.min, id="Treatment")  # convert to long format
rm(datpOpt.feature.min)
# total counts and indexed counts per cpd [9 reps per cmpd for every feature]
datpOpt.feature.min.long$id=1; datpOpt.feature.min.long<-data.table(datpOpt.feature.min.long)

datpOpt.feature.min.long[, totalreplicates := .N, by = c('Treatment','variable')]
datpOpt.feature.min.long[, indexreplicates := cumsum(id), by = c('Treatment','variable')]
datpOpt.feature.min.long=datpOpt.feature.min.long[,-c('id','totalreplicates')]

#renaming
datpOpt.feature.min.long=
  plyr::rename(datpOpt.feature.min.long,c('Treatment'='Treatment','variable'='feature','value'='value'))

# rename TRT levels
datpOpt.feature.min.long$Treatment=  gsub("in Colon-CytoSkeleton-Endosomes","",datpOpt.feature.min.long$Treatment)

#plots of replicates
datpOpt.feature.min.long$Treatment=as.factor(datpOpt.feature.min.long$Treatment)
replicatesall_featselectedmin=  ggplot(data = datpOpt.feature.min.long,aes(x = feature, y = value,group=interaction(indexreplicates,Treatment),col=Treatment)) +
  geom_point(pch=19, cex=0.55) + theme_bw() + geom_line(lwd=0.6)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab('feature')+
  ylab('replicate value')+ guides(fill=guide_legend(ncol=2)) + theme(legend.position="top",legend.text=element_text(size=10)) 

ggsave('report/dataView/EDA/EDA_minimal_diffTRT.pdf',replicatesall_featselectedmin)

rm(datpOpt.feature.min.long,datpOpt.feature.max.long)

### DATASET2
eda2<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/dataView/EDA'

# 650 treatments, optimal features [untransformed data]
load('data/p1602xxPPS_means_preprocessed.Rdata')
datpOpt.feature<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features_mRMR_optN"))]
rm(list=c('datp'))

# plot profiles for treatments that are highly similar via Hotelling's T2 statistic+those that are very dissimilar using T2
# idea from the Lwafula_HotellingT2_calculation-Final script
# [cmpd892 @ 9 uM in Hepa-ER-Mito & cmpd94 @ 9 uM in Hepa-ER-Mito with T2=242201.3] ; tofauti kabisa
# [cmpd417 @ 1 uM in Hepa-ER-Mito & cmpd417 @ 3 uM in Hepa-ER-Mito with T2=2.42051] ; tofauti kidogo

trtkp<-c('cmpd892 @ 9 uM in Hepa-ER-Mito','cmpd94 @ 9 uM in Hepa-ER-Mito','cmpd417 @ 1 uM in Hepa-ER-Mito','cmpd417 @ 3 uM in Hepa-ER-Mito')

# keep only the target ones
datpOpt.feature<-datpOpt.feature[datpOpt.feature[['Treatment']] %in% trtkp,]

# maximal difference combination
datpOpt.feature.max<-
  datpOpt.feature[datpOpt.feature[['Treatment']] %in% c('cmpd892 @ 9 uM in Hepa-ER-Mito','cmpd94 @ 9 uM in Hepa-ER-Mito'),]

datpOpt.feature.max.long <- melt(datpOpt.feature.max, id="Treatment")  # convert to long format
rm(datpOpt.feature.max)
# total counts and indexed counts per cpd [9 reps per cmpd for every feature]
datpOpt.feature.max.long$id=1; datpOpt.feature.max.long<-data.table(datpOpt.feature.max.long)

datpOpt.feature.max.long[, totalreplicates := .N, by = c('Treatment','variable')]
datpOpt.feature.max.long[, indexreplicates := cumsum(id), by = c('Treatment','variable')]
datpOpt.feature.max.long=datpOpt.feature.max.long[,-c('id','totalreplicates')]

#renaming
datpOpt.feature.max.long=
  plyr::rename(datpOpt.feature.max.long,c('Treatment'='Treatment','variable'='feature','value'='value'))

# rename TRT levels
datpOpt.feature.max.long$Treatment=  gsub("in Hepa-ER-Mito","",datpOpt.feature.max.long$Treatment)

#plots of replicates
datpOpt.feature.max.long$Treatment=as.factor(datpOpt.feature.max.long$Treatment)
replicatesall_featselected=  ggplot(data = datpOpt.feature.max.long,aes(x = feature, y = value,group=interaction(indexreplicates,Treatment),col=Treatment)) +
  geom_point(pch=19, cex=0.55) + theme_bw() + geom_line(lwd=0.6)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab('feature')+
  ylab('replicate value')+ guides(fill=guide_legend(ncol=2)) + theme(legend.position="top",legend.text=element_text(size=10)) 

ggsave('E:/uhasselt/22yrsem/MThesis/DATASET2/report/dataView/EDA/D2_EDA_maximal_diffTRT.pdf',replicatesall_featselected)

# minimal diff
datpOpt.feature.min<-
  datpOpt.feature[datpOpt.feature[['Treatment']] %in% c('cmpd417 @ 1 uM in Hepa-ER-Mito','cmpd417 @ 3 uM in Hepa-ER-Mito'),]

datpOpt.feature.min.long <- melt(datpOpt.feature.min, id="Treatment")  # convert to long format
rm(datpOpt.feature.min)
# total counts and indexed counts per cpd [9 reps per cmpd for every feature]
datpOpt.feature.min.long$id=1; datpOpt.feature.min.long<-data.table(datpOpt.feature.min.long)

datpOpt.feature.min.long[, totalreplicates := .N, by = c('Treatment','variable')]
datpOpt.feature.min.long[, indexreplicates := cumsum(id), by = c('Treatment','variable')]
datpOpt.feature.min.long=datpOpt.feature.min.long[,-c('id','totalreplicates')]

#renaming
datpOpt.feature.min.long=
  plyr::rename(datpOpt.feature.min.long,c('Treatment'='Treatment','variable'='feature','value'='value'))

# rename TRT levels
datpOpt.feature.min.long$Treatment=  gsub("in Hepa-ER-Mito","",datpOpt.feature.min.long$Treatment)

#plots of replicates
datpOpt.feature.min.long$Treatment=as.factor(datpOpt.feature.min.long$Treatment)
replicatesall_featselectedmin=  ggplot(data = datpOpt.feature.min.long,aes(x = feature, y = value,group=interaction(indexreplicates,Treatment),col=Treatment)) +
  geom_point(pch=19, cex=0.55) + theme_bw() + geom_line(lwd=0.6)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab('feature')+
  ylab('replicate value')+ guides(fill=guide_legend(ncol=2)) + theme(legend.position="top",legend.text=element_text(size=10)) 

ggsave('DATASET2/report/dataView/EDA/D2_EDA_minimal_diffTRT.pdf',replicatesall_featselectedmin)



# distribution of features b4 and after transformation [only features selected, active treatments]
load('data/p1601xxPPS_means_preprocessed.Rdata')
datpOpt.feature<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features_mRMR_optN"))]

features_opts_orig<-attr(datp, "features_mRMR_optN"); features_opts <- gsub(".mean.zGSLC", "",features_opts_orig); features_opts
#remove the .mean.zGSLC
names(datpOpt.feature) <- gsub(".mean.zGSLC", "", names(datpOpt.feature))
rm(list=c('datp'))

# histograms

myhist <- hist(datpOpt.feature[,features_opts[10]],plot = FALSE)
multiplier <- myhist$counts/ myhist$density
mydensity <- density(datpOpt.feature[,features_opts[10]])
mydensity$y <- mydensity$y * multiplier[1]
mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
mydensity$gr<-'untransformed'


plot(mydensity[mydensity$gr=='untransformed',]$y~mydensity[mydensity$gr=='untransformed',]$x,type='l',lwd=1,lty=2,ylab='Density',
     main=eval(bquote(expression(paste('Distribution of'~.(features_opts[10])~'b4 transformation')))), xlab='value')

#histograms for all features selected in the untransformed cases
pdf('report/dataView/EDA/data1_distr_1_9.pdf')
par(mfrow=c(3,3))
for(i in 1:length(features_opts[1:9])){
myhist <- hist(datpOpt.feature[,features_opts[i]],plot = FALSE)
multiplier <- myhist$counts / myhist$density
mydensity <- density(datpOpt.feature[,features_opts[i]])
mydensity$y <- mydensity$y * multiplier[1]
mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
mydensity$gr<-'untransformed'

plot(mydensity$y~mydensity$x,type='l',lwd=1,lty=1,ylab='Density',xlab='value')
title(eval(bquote(expression(paste('Distribution of'~.(features_opts[i]))))), sub = " ",
      cex.main = 1,   font.main=2, col.main= "black",
      cex.sub = 0.75, font.sub = 3, col.sub = "black")
}
dev.off()

#histograms for all features selected in the untransformed cases
pdf('report/dataView/EDA/data1_distr_10_18.pdf')
par(mfrow=c(3,3))
for(i in 1:length(features_opts[10:18])){
  myhist <- hist(datpOpt.feature[,features_opts[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,features_opts[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  plot(mydensity$y~mydensity$x,type='l',lwd=1,lty=1,ylab='Density',xlab='value')
  title(eval(bquote(expression(paste('Distribution of'~.(features_opts[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
}
dev.off()

pdf('report/dataView/EDA/data1_distr_19.pdf')
par(mfrow=c(1,1))
  myhist <- hist(datpOpt.feature[,features_opts[19]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,features_opts[19]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  plot(mydensity$y~mydensity$x,type='l',lwd=1,lty=1,ylab='Density',xlab='value')
  title(eval(bquote(expression(paste('Distribution of'~.(features_opts[19]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
dev.off()

############

# transformed data for same features for lambda=10.5
load('data/glog transformed datasets by lambda - Copy/p1601xxPPS_means_preprocessedfeat_act_10.5.Rdata')
glogdatpOpt.feature<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment",  attr(glogdatp, "features_mRMR_optN"))]
features_glogopts<-attr(glogdatp, "features_mRMR_optN")
features_glogopts <- gsub(".zGSLC", "",features_glogopts); features_glogopts
#remove the .mean.zGSLC
names(glogdatpOpt.feature) <- gsub(".zGSLC", "", names(glogdatpOpt.feature))

rm(list=c('glogdatp'))

# distributions
#histograms for all features selected in the transformed cases
pdf('report/dataView/EDA/data1_transformed_distr_1_12.pdf')
par(mfrow=c(4,3))
for(i in 1:length(features_glogopts[1:12])){
  myhist <- hist(glogdatpOpt.feature[,features_glogopts[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(glogdatpOpt.feature[,features_glogopts[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  
  plot(mydensity$y~mydensity$x,type='l',lwd=1,lty=1,ylab='Density',xlab='value')
  title(eval(bquote(expression(paste('Distribution of'~.(features_glogopts[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
}
dev.off()

pdf('report/dataView/EDA/data1_transformed_distr_13_21.pdf')
par(mfrow=c(3,3))
for(i in 1:length(features_opts[13:21])){
  myhist <- hist(glogdatpOpt.feature[,features_glogopts[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(glogdatpOpt.feature[,features_glogopts[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  
  plot(mydensity$y~mydensity$x,type='l',lwd=1,lty=1,ylab='Density',xlab='value')
  title(eval(bquote(expression(paste('Distribution of'~.(features_glogopts[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")}
dev.off()

# distribution of features b4 and after transformation [those that are common in both data sets]

rm(list = ls())
datawd<-'E:/uhasselt/22yrsem/MThesis/data'
eda<-'E:/uhasselt/22yrsem/MThesis/report/dataView/EDA'

# untransformed
load('E:/uhasselt/22yrsem/MThesis/data/p1601xxPPS_means_preprocessed.Rdata')
datpOpt.feature<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features_mRMR_optN"))]

features_opts_orig<-attr(datp, "features_mRMR_optN"); features_opts <- gsub(".mean.zGSLC", "",features_opts_orig); features_opts
#remove the .mean.zGSLC
names(datpOpt.feature) <- gsub(".mean.zGSLC", "", names(datpOpt.feature))
rm(list=c('datp'))

#transformed, lambda= 10.5
load('data/glog transformed datasets by lambda/p1601xxPPS_means_preprocessedfeat_act_10.5.Rdata')
glogdatpOpt.feature<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(glogdatp, "features_mRMR_optN"))]

glogfeatures_opts_orig<-attr(glogdatp, "features_mRMR_optN")
glogfeatures_opts <- gsub(".zGSLC", "",glogfeatures_opts_orig); glogfeatures_opts
#remove the .mean.zGSLC
names(glogdatpOpt.feature) <- gsub(".zGSLC", "", names(glogdatpOpt.feature))
rm(list=c('glogdatp'))

# columns that are both in untransformed and transformed
glogfeatureint<-intersect(features_opts,glogfeatures_opts)
# save
pdf<-'report/dataView/EDA/data1_distr_overlaps_pre_post_transformation.pdf'
pdf(pdf)
par(mfrow=c(2,2))
for (i in 1:length(glogfeatureint)){
  myhist <- hist(datpOpt.feature[,glogfeatureint[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,glogfeatureint[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  myhist0 <- hist(glogdatpOpt.feature[,glogfeatureint[i]],plot = FALSE)
  multiplier <- myhist0$counts / myhist0$density
  mydensity0 <- density(glogdatpOpt.feature[,glogfeatureint[i]])
  mydensity0$y <- mydensity0$y * multiplier[1]
  mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
  mydensity0$gr<-'transformed'
  mydensity<-data.frame(rbind(mydensity,mydensity0))
  
  par(mar = c(5,5,2,5))
  with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                      xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
  par(new = T)
  with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                      axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
  axis(side = 4,col='blue')
  mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
  title(eval(bquote(expression(paste('Distribution of '~.(glogfeatureint[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
  
}
dev.off()

### FEATURES ABSENT IN SELECTED FEATURES PRE-TRANSFORMATION, BUT SELECTED POST-TRANSFORMATION: CHECK THAT OF 1ST FEATURE IN 9_12.PDF


# `````````````` OBJ.1 : showing the differentiating ability of the selected features ``````````````````````````````````````````
# prove that the glog transformation did little to improve distributions of features pe- and post-transformation [using features not 
# selected pre-transformation but selected after]
# LWafula 27072017
# Thz script is primarily for all EDA for all data-performed on all the cellLines
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
library(Hotelling) # for Hotelling calculation
library(qpcR) # cbinding unequal lengths
# distribution of features not selected b4 but selected after transformation [those in transformed but not in transformed]

rm(list = ls())
datawd<-'E:/uhasselt/22yrsem/MThesis/data'
eda<-'E:/uhasselt/22yrsem/MThesis/report/dataView/EDA'

# untransformed
load('E:/uhasselt/22yrsem/MThesis/data/p1601xxPPS_means_preprocessed.Rdata')
datpOpt.feature<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features_mRMR_optN"))]

features_opts_orig<-attr(datp, "features_mRMR_optN"); features_opts <- gsub(".mean.zGSLC", "",features_opts_orig); features_opts
#remove the .mean.zGSLC
names(datpOpt.feature) <- gsub(".mean.zGSLC", "", names(datpOpt.feature))

#transformed, lambda= 10.5
load('data/glog transformed datasets by lambda/p1601xxPPS_means_preprocessedfeat_act_10.5.Rdata')
glogdatpOpt.feature<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(glogdatp, "features_mRMR_optN"))]

glogfeatures_opts_orig<-attr(glogdatp, "features_mRMR_optN")
glogfeatures_opts <- gsub(".zGSLC", "",glogfeatures_opts_orig); glogfeatures_opts
#remove the .mean.zGSLC
names(glogdatpOpt.feature) <- gsub(".zGSLC", "", names(glogdatpOpt.feature))
rm(list=c('glogdatp'))

# columns that are in transformed and not in untransformed
# ftes <- qpcR:::cbind.na(sort(glogfeatures_opts),sort(features_opts))
glogfeaturediff<-setdiff(glogfeatures_opts,features_opts)
# save

# for both datasets, keep only those in transformed but not in untransformed
glogdatpOpt.feature<-glogdatpOpt.feature[,colnames(glogdatpOpt.feature) %in% glogfeaturediff]
datpOpt.feature<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features"))]
names(datpOpt.feature) <- gsub(".mean.zGSLC", "", names(datpOpt.feature))
datpOpt.feature<-datpOpt.feature[,colnames(datpOpt.feature) %in% glogfeaturediff]
rm(list=c('datp'))

pdf<-'report/dataView/EDA/data1_distr_setdiff_pre_post_transformation_1to4.pdf'
pdf(pdf)
par(mfrow=c(2,2))
for (i in 1:4){
  myhist <- hist(datpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,glogfeaturediff[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  myhist0 <- hist(glogdatpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist0$counts / myhist0$density
  mydensity0 <- density(glogdatpOpt.feature[,glogfeaturediff[i]])
  mydensity0$y <- mydensity0$y * multiplier[1]
  mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
  mydensity0$gr<-'transformed'
  mydensity<-data.frame(rbind(mydensity,mydensity0))
  
  par(mar = c(5,5,2,5))
  with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                      xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
  par(new = T)
  with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                      axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
  axis(side = 4,col='blue')
  mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
  title(eval(bquote(expression(paste('Distribution of '~.(glogfeaturediff[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
  
}
dev.off()

pdf<-'report/dataView/EDA/data1_distr_setdiff_pre_post_transformation_5to8.pdf'
pdf(pdf)
par(mfrow=c(2,2))
for (i in 5:8){
  myhist <- hist(datpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,glogfeaturediff[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  myhist0 <- hist(glogdatpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist0$counts / myhist0$density
  mydensity0 <- density(glogdatpOpt.feature[,glogfeaturediff[i]])
  mydensity0$y <- mydensity0$y * multiplier[1]
  mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
  mydensity0$gr<-'transformed'
  mydensity<-data.frame(rbind(mydensity,mydensity0))
  
  par(mar = c(5,5,2,5))
  with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                      xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
  par(new = T)
  with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                      axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
  axis(side = 4,col='blue')
  mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
  title(eval(bquote(expression(paste('Distribution of '~.(glogfeaturediff[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
  
}
dev.off()

pdf<-'report/dataView/EDA/data1_distr_setdiff_pre_post_transformation_9to12.pdf'
pdf(pdf)
par(mfrow=c(2,2))
for (i in 9:12){
  myhist <- hist(datpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,glogfeaturediff[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  myhist0 <- hist(glogdatpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist0$counts / myhist0$density
  mydensity0 <- density(glogdatpOpt.feature[,glogfeaturediff[i]])
  mydensity0$y <- mydensity0$y * multiplier[1]
  mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
  mydensity0$gr<-'transformed'
  mydensity<-data.frame(rbind(mydensity,mydensity0))
  
  par(mar = c(5,5,2,5))
  with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                      xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
  par(new = T)
  with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                      axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
  axis(side = 4,col='blue')
  mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
  title(eval(bquote(expression(paste('Distribution of '~.(glogfeaturediff[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
  
}
dev.off()

pdf<-'report/dataView/EDA/data1_distr_setdiff_pre_post_transformation_13to16.pdf'
pdf(pdf)
par(mfrow=c(2,2))
for (i in 13:16){
  myhist <- hist(datpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,glogfeaturediff[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  myhist0 <- hist(glogdatpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist0$counts / myhist0$density
  mydensity0 <- density(glogdatpOpt.feature[,glogfeaturediff[i]])
  mydensity0$y <- mydensity0$y * multiplier[1]
  mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
  mydensity0$gr<-'transformed'
  mydensity<-data.frame(rbind(mydensity,mydensity0))
  
  par(mar = c(5,5,2,5))
  with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                      xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
  par(new = T)
  with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                      axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
  axis(side = 4,col='blue')
  mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
  title(eval(bquote(expression(paste('Distribution of '~.(glogfeaturediff[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
  
}
dev.off()

# 17th feature
pdf<-'report/dataView/EDA/data1_distr_setdiff_pre_post_transformation_17.pdf'
pdf(pdf)
par(mfrow=c(1,1))
myhist <- hist(datpOpt.feature[,glogfeaturediff[17]],plot = FALSE)
multiplier <- myhist$counts / myhist$density
mydensity <- density(datpOpt.feature[,glogfeaturediff[17]])
mydensity$y <- mydensity$y * multiplier[1]
mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
mydensity$gr<-'untransformed'

myhist0 <- hist(glogdatpOpt.feature[,glogfeaturediff[17]],plot = FALSE)
multiplier <- myhist0$counts / myhist0$density
mydensity0 <- density(glogdatpOpt.feature[,glogfeaturediff[17]])
mydensity0$y <- mydensity0$y * multiplier[1]
mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
mydensity0$gr<-'transformed'
mydensity<-data.frame(rbind(mydensity,mydensity0))

par(mar = c(5,5,2,5))
with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                    xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
par(new = T)
with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                    axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
axis(side = 4,col='blue')
mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
title(eval(bquote(expression(paste('Distribution of '~.(glogfeaturediff[i]))))), sub = " ",
      cex.main = 1,   font.main=2, col.main= "black",
      cex.sub = 0.75, font.sub = 3, col.sub = "black")
dev.off()
##################

# FEATURES COMMON IN BOTH UNTRANSFORMED AND TRANSFORMED DATA, LAMBDA= 0.5
# features in DATASET 2
rm(list = ls())
datawd<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data'
eda<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/dataView/EDA'

# untransformed
load('E:/uhasselt/22yrsem/MThesis/data/p1602xxPPS_means_preprocessed.Rdata')
datpOpt.feature<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features_mRMR_optN"))]

features_opts_orig<-attr(datp, "features_mRMR_optN"); features_opts <- gsub(".mean.zGSLC", "",features_opts_orig); features_opts
#remove the .mean.zGSLC
names(datpOpt.feature) <- gsub(".mean.zGSLC", "", names(datpOpt.feature))
rm(list=c('datp'))

#transformed, lambda= 0.5
load('DATASET2/data/glog transformed datasets by lambda/p1602xxPPS_means_preprocessedfeat_act_0.5.Rdata')
glogdatpOpt.feature<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(glogdatp, "features_mRMR_optN"))]

glogfeatures_opts_orig<-attr(glogdatp, "features_mRMR_optN")
glogfeatures_opts <- gsub(".zGSLC", "",glogfeatures_opts_orig); glogfeatures_opts
#remove the .mean.zGSLC
names(glogdatpOpt.feature) <- gsub(".zGSLC", "", names(glogdatpOpt.feature))
rm(list=c('glogdatp'))

# columns that are both in untransformed and transformed
glogfeatureint<-intersect(features_opts,glogfeatures_opts)

# save
pdf<-'report/dataView/EDA/data2_distr_overlaps_pre_post_transformation1_4.pdf'
pdf(pdf)
par(mfrow=c(2,2))
for (i in 1:4){
  myhist <- hist(datpOpt.feature[,glogfeatureint[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,glogfeatureint[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  myhist0 <- hist(glogdatpOpt.feature[,glogfeatureint[i]],plot = FALSE)
  multiplier <- myhist0$counts / myhist0$density
  mydensity0 <- density(glogdatpOpt.feature[,glogfeatureint[i]])
  mydensity0$y <- mydensity0$y * multiplier[1]
  mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
  mydensity0$gr<-'transformed'
  mydensity<-data.frame(rbind(mydensity,mydensity0))
  
  par(mar = c(5,5,2,5))
  with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                      xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
  par(new = T)
  with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                      axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
  axis(side = 4,col='blue')
  mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
  title(eval(bquote(expression(paste('Distribution of '~.(glogfeatureint[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
  
}
dev.off()

# save
pdf<-'report/dataView/EDA/data2_distr_overlaps_pre_post_transformation5_8.pdf'
pdf(pdf)
par(mfrow=c(2,2))
for (i in 5:8){
  myhist <- hist(datpOpt.feature[,glogfeatureint[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,glogfeatureint[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  myhist0 <- hist(glogdatpOpt.feature[,glogfeatureint[i]],plot = FALSE)
  multiplier <- myhist0$counts / myhist0$density
  mydensity0 <- density(glogdatpOpt.feature[,glogfeatureint[i]])
  mydensity0$y <- mydensity0$y * multiplier[1]
  mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
  mydensity0$gr<-'transformed'
  mydensity<-data.frame(rbind(mydensity,mydensity0))
  
  par(mar = c(5,5,2,5))
  with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                      xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
  par(new = T)
  with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                      axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
  axis(side = 4,col='blue')
  mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
  title(eval(bquote(expression(paste('Distribution of '~.(glogfeatureint[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
  
}
dev.off()


# save
pdf<-'report/dataView/EDA/data2_distr_overlaps_pre_post_transformation9_12.pdf'
pdf(pdf)
par(mfrow=c(2,2))
for (i in 9:12){
  myhist <- hist(datpOpt.feature[,glogfeatureint[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,glogfeatureint[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  myhist0 <- hist(glogdatpOpt.feature[,glogfeatureint[i]],plot = FALSE)
  multiplier <- myhist0$counts / myhist0$density
  mydensity0 <- density(glogdatpOpt.feature[,glogfeatureint[i]])
  mydensity0$y <- mydensity0$y * multiplier[1]
  mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
  mydensity0$gr<-'transformed'
  mydensity<-data.frame(rbind(mydensity,mydensity0))
  
  par(mar = c(5,5,2,5))
  with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                      xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
  par(new = T)
  with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                      axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
  axis(side = 4,col='blue')
  mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
  title(eval(bquote(expression(paste('Distribution of '~.(glogfeatureint[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
  
}
dev.off()

pdf<-'report/dataView/EDA/data2_distr_overlaps_pre_post_transformation13_16.pdf'
pdf(pdf)
par(mfrow=c(2,2))
for (i in 13:16){
  myhist <- hist(datpOpt.feature[,glogfeatureint[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,glogfeatureint[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  myhist0 <- hist(glogdatpOpt.feature[,glogfeatureint[i]],plot = FALSE)
  multiplier <- myhist0$counts / myhist0$density
  mydensity0 <- density(glogdatpOpt.feature[,glogfeatureint[i]])
  mydensity0$y <- mydensity0$y * multiplier[1]
  mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
  mydensity0$gr<-'transformed'
  mydensity<-data.frame(rbind(mydensity,mydensity0))
  
  par(mar = c(5,5,2,5))
  with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                      xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
  par(new = T)
  with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                      axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
  axis(side = 4,col='blue')
  mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
  title(eval(bquote(expression(paste('Distribution of '~.(glogfeatureint[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
  
}
dev.off()



### FEATURES ABSENT IN SELECTED FEATURES PRE-TRANSFORMATION, BUT SELECTED POST-TRANSFORMATION: CHECK THAT OF 1ST FEATURE IN 
## 1,3 , 6,7,8, and 18: right and in  some cases left skews present pre-transformation has been in most cases minimized, not that much 
# but yes, something happened


# `````````````` OBJ.1 : showing the differentiating ability of the selected features ``````````````````````````````````````````
# prove that the glog transformation did little to improve distributions of features pe- and post-transformation [using features not 
# selected pre-transformation but selected after]
# LWafula 27072017
# Thz script is primarily for all EDA for all data-performed on all the cellLines
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
library(Hotelling) # for Hotelling calculation
library(qpcR) # cbinding unequal lengths
# distribution of features not selected b4 but selected after transformation [those in transformed but not in transformed]

rm(list = ls())
datawd<-'E:/uhasselt/22yrsem/MThesis/data'
eda<-'E:/uhasselt/22yrsem/MThesis/report/dataView/EDA'

# untransformed
load('E:/uhasselt/22yrsem/MThesis/data/p1602xxPPS_means_preprocessed.Rdata')
datpOpt.feature<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features_mRMR_optN"))]

features_opts_orig<-attr(datp, "features_mRMR_optN"); features_opts <- gsub(".mean.zGSLC", "",features_opts_orig); features_opts
#remove the .mean.zGSLC
names(datpOpt.feature) <- gsub(".mean.zGSLC", "", names(datpOpt.feature))

#transformed, lambda= 0.5
load('data/glog transformed datasets by lambda/p1602xxPPS_means_preprocessedfeat_act_0.5.Rdata')
glogdatpOpt.feature<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(glogdatp, "features_mRMR_optN"))]

glogfeatures_opts_orig<-attr(glogdatp, "features_mRMR_optN")
glogfeatures_opts <- gsub(".zGSLC", "",glogfeatures_opts_orig); glogfeatures_opts
#remove the .mean.zGSLC
names(glogdatpOpt.feature) <- gsub(".zGSLC", "", names(glogdatpOpt.feature))
rm(list=c('glogdatp'))

# columns that are in transformed and not in untransformed
# ftes <- qpcR:::cbind.na(sort(glogfeatures_opts),sort(features_opts))
glogfeaturediff<-setdiff(glogfeatures_opts,features_opts)
# save

# for both datasets, keep only those in transformed but not in untransformed
glogdatpOpt.feature<-glogdatpOpt.feature[,colnames(glogdatpOpt.feature) %in% glogfeaturediff]
datpOpt.feature<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features"))]
names(datpOpt.feature) <- gsub(".mean.zGSLC", "", names(datpOpt.feature))
datpOpt.feature<-datpOpt.feature[,colnames(datpOpt.feature) %in% glogfeaturediff]
rm(list=c('datp'))

pdf<-'report/dataView/EDA/data2_distr_setdiff_pre_post_transformation_1to4.pdf'
pdf(pdf)
par(mfrow=c(2,2))
for (i in 1:4){
  myhist <- hist(datpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,glogfeaturediff[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  myhist0 <- hist(glogdatpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist0$counts / myhist0$density
  mydensity0 <- density(glogdatpOpt.feature[,glogfeaturediff[i]])
  mydensity0$y <- mydensity0$y * multiplier[1]
  mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
  mydensity0$gr<-'transformed'
  mydensity<-data.frame(rbind(mydensity,mydensity0))
  
  par(mar = c(5,5,2,5))
  with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                      xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
  par(new = T)
  with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                      axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
  axis(side = 4,col='blue')
  mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
  title(eval(bquote(expression(paste('Distribution of '~.(glogfeaturediff[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
  
}
dev.off()

pdf<-'report/dataView/EDA/data2_distr_setdiff_pre_post_transformation_5to8.pdf'
pdf(pdf)
par(mfrow=c(2,2))
for (i in 5:8){
  myhist <- hist(datpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,glogfeaturediff[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  myhist0 <- hist(glogdatpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist0$counts / myhist0$density
  mydensity0 <- density(glogdatpOpt.feature[,glogfeaturediff[i]])
  mydensity0$y <- mydensity0$y * multiplier[1]
  mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
  mydensity0$gr<-'transformed'
  mydensity<-data.frame(rbind(mydensity,mydensity0))
  
  par(mar = c(5,5,2,5))
  with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                      xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
  par(new = T)
  with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                      axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
  axis(side = 4,col='blue')
  mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
  title(eval(bquote(expression(paste('Distribution of '~.(glogfeaturediff[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
  
}
dev.off()

pdf<-'report/dataView/EDA/data2_distr_setdiff_pre_post_transformation_9to12.pdf'
pdf(pdf)
par(mfrow=c(2,2))
for (i in 9:12){
  myhist <- hist(datpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,glogfeaturediff[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  myhist0 <- hist(glogdatpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist0$counts / myhist0$density
  mydensity0 <- density(glogdatpOpt.feature[,glogfeaturediff[i]])
  mydensity0$y <- mydensity0$y * multiplier[1]
  mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
  mydensity0$gr<-'transformed'
  mydensity<-data.frame(rbind(mydensity,mydensity0))
  
  par(mar = c(5,5,2,5))
  with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                      xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
  par(new = T)
  with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                      axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
  axis(side = 4,col='blue')
  mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
  title(eval(bquote(expression(paste('Distribution of '~.(glogfeaturediff[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
  
}
dev.off()

pdf<-'report/dataView/EDA/data2_distr_setdiff_pre_post_transformation_13to16.pdf'
pdf(pdf)
par(mfrow=c(2,2))
for (i in 13:16){
  myhist <- hist(datpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,glogfeaturediff[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  myhist0 <- hist(glogdatpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist0$counts / myhist0$density
  mydensity0 <- density(glogdatpOpt.feature[,glogfeaturediff[i]])
  mydensity0$y <- mydensity0$y * multiplier[1]
  mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
  mydensity0$gr<-'transformed'
  mydensity<-data.frame(rbind(mydensity,mydensity0))
  
  par(mar = c(5,5,2,5))
  with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                      xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
  par(new = T)
  with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                      axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
  axis(side = 4,col='blue')
  mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
  title(eval(bquote(expression(paste('Distribution of '~.(glogfeaturediff[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
  
}
dev.off()

# 17th and 18th features
pdf<-'report/dataView/EDA/data2_distr_setdiff_pre_post_transformation_17_18.pdf'
pdf(pdf)
par(mfrow=c(2,2))
for (i in 17:18){
  myhist <- hist(datpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(datpOpt.feature[,glogfeaturediff[i]])
  mydensity$y <- mydensity$y * multiplier[1]
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'untransformed'
  
  myhist0 <- hist(glogdatpOpt.feature[,glogfeaturediff[i]],plot = FALSE)
  multiplier <- myhist0$counts / myhist0$density
  mydensity0 <- density(glogdatpOpt.feature[,glogfeaturediff[i]])
  mydensity0$y <- mydensity0$y * multiplier[1]
  mydensity0<-data.frame(cbind(y=mydensity0$y,x=mydensity0$x))
  mydensity0$gr<-'transformed'
  mydensity<-data.frame(rbind(mydensity,mydensity0))
  
  par(mar = c(5,5,2,5))
  with(mydensity,plot(mydensity[mydensity$gr=='untransformed',]$x,mydensity[mydensity$gr=='untransformed',]$y,ylab='untransformed',
                      xlab='value',lty=1,lwd=2,cex=0.5,type='l'))
  par(new = T)
  with(mydensity,plot(mydensity[mydensity$gr=='transformed',]$x,mydensity[mydensity$gr=='transformed',]$y, 
                      axes=F, xlab=NA, ylab=NA, cex=0.6,type='l',col='blue',lwd=2))
  axis(side = 4,col='blue')
  mtext(side = 4, line = 3,eval(bquote(expression(paste('transformed')))),,col='blue')
  title(eval(bquote(expression(paste('Distribution of '~.(glogfeaturediff[i]))))), sub = " ",
        cex.main = 1,   font.main=2, col.main= "black",
        cex.sub = 0.75, font.sub = 3, col.sub = "black")
}
dev.off()



###########################################################################

############ are all active treatments in untransformed picked in those transformations that seem to perform well? NOP
# data set 1
rm(list=ls())
load('data/p1601xxPPS_means_preprocessed.Rdata')
datpOpt.feature<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features_mRMR_optN"))]
rm(list=c('datp'))
dim(table(datpOpt.feature$Treatment))

# treatments b4 transformations
trtsb4<-names(table(datpOpt.feature$Treatment))

## treatments after transformation (lambda= 10.5)
load('data/glog transformed datasets by lambda/p1601xxPPS_means_preprocessedfeat_act_10.5.RData')
glogdatpOpt.feature<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(glogdatp, "features_mRMR_optN"))]
rm(list=c('glogdatp'))
dim(table(glogdatpOpt.feature$Treatment))
trts.after105<-names(table(glogdatpOpt.feature$Treatment))

# set difference not equal to dim(table(glogdatpOpt.feature$Treatment))- dim(table(datpOpt.feature$Treatment))
c(dim(table(glogdatpOpt.feature$Treatment)), dim(table(datpOpt.feature$Treatment)),
  dim(table(glogdatpOpt.feature$Treatment))- dim(table(datpOpt.feature$Treatment)),length(setdiff(trts.after105,trtsb4)))

setdiff(trtsb4,trts.after105) # 24 in untransformed that don't make it in transformed
setdiff(trts.after105,trtsb4) # 20 new ones in tranformed compared to untransformed

# common in both pre-and post-transformation
intersect(trtsb4,trts.after105) # 683

## data set 2
load('data/p1602xxPPS_means_preprocessed.Rdata')
dat2pOpt.feature<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features_mRMR_optN"))]
rm(list=c('datp'))
dim(table(dat2pOpt.feature$Treatment))

# treatments b4 transformations
trtsb4.d2<-names(table(dat2pOpt.feature$Treatment))

## treatments after transformation (lambda= 0.5)
load('DATASET2/data/glog transformed datasets by lambda/p1602xxPPS_means_preprocessedfeat_act_0.5.RData')
glogdat2pOpt.feature<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(glogdatp, "features_mRMR_optN"))]
rm(list=c('glogdatp'))
dim(table(glogdat2pOpt.feature$Treatment))
trts.after05<-names(table(glogdat2pOpt.feature$Treatment))

# set difference not equal to dim(table(glogdatpOpt.feature$Treatment))- dim(table(datpOpt.feature$Treatment))
c(dim(table(glogdat2pOpt.feature$Treatment)), dim(table(dat2pOpt.feature$Treatment)),
  dim(table(glogdat2pOpt.feature$Treatment))- dim(table(dat2pOpt.feature$Treatment)),length(setdiff(trts.after05,trtsb4.d2)))

setdiff(trtsb4.d2,trts.after05) # 8 in untransformed that don't make it in transformed
setdiff(trts.after05,trtsb4.d2) # 20 new ones in transformed compared to untransformed

# common actively-called treatments
intersect(trts.after05,trtsb4.d2) #642














