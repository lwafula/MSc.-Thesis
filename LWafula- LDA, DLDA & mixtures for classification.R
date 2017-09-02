# LDA/ DLDA on the untransformed data

# start with two treatments  
#
# LWafula : 07-06-2017
#Run- the data preprocessing script [doesn't take long]

# for the gereral lambda to check the actively-called per lambda
# and also the number of features selected
# data processing from single-cell data to well-level data

# data viewing
# LWafula: 02/06/2017
rm(list=ls())

setwd('E:/uhasselt/22yrsem/MThesis')
library(hexbin)
library(ggplot2)
library(matlab)
library(plyr)
library(dplyr)
library(data.table)
library(flux)
library(WGCNA) #https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/InstallationInstructions.html
library(knitr) #for kable
library("reshape2")
library(MASS)
# plots for LDA
# https://gist.github.com/thigm85/8424654
require(ggplot2)
require(scales)
require(gridExtra)

# data wd
datawd<-'E:/uhasselt/22yrsem/MThesis/data'

# untransformed feature selected dataset
#keep the selected features only, and select 2 active treatments preferably from different compounds so that we see the discrimination
datpfeatact<-'p1601xxPPS_means_preprocessedfeat_act.Rdata'
datpfeatact<-file.path(datawd,datpfeatact)
load(datpfeatact)

#remove the endosome thing in the treatment name
datpfeatact$Treatment=   gsub("in Colon-CytoSkeleton-Endosomes","",datpfeatact$Treatment)

#candidate compounds sent by Steffen
cmpd<-c('cmpd291', 'cmpd788', 'cmpd653', 'cmpd626', 'cmpd104', 'cmpd262', 'cmpd298', 'cmpd595', 'cmpd464',
        'cmpd275', 'cmpd306', 'cmpd306', 'cmpd486', 'cmpd486')
# as sent
# treatset<-c('cmpd291 @ 3 uM','cmpd788 @ 3 uM','cmpd653 @ 0.334 uM',
#             'cmpd626 @ 3 uM','cmpd104 @ 3 uM','cmpd262 @ 3 uM','cmpd298 @ 3 uM',
#             'cmpd595 @ 1 uM','cmpd464 @ 9 uM','cmpd275 @ 1 uM','cmpd306 @ 1 uM',
#             'cmpd306 @ 3 uM','cmpd486 @ 1 uM','cmpd486 @ 3 uM')

#keep relevant cmpds and concentration levels
datpfeatact9<-datpfeatact[CMPD_ID %in% cmpd,]

#treats to keep 
treatkp<-names(table(datpfeatact9$Treatment)[c(11,42,36, 34,1,3,15,30,24,6,18,19,26,27)])
datpfeatact9<-datpfeatact[Treatment %in% treatkp,]

#plot two features and separate by Treatment [2]
datpfeatact9<-data.frame(datpfeatact9)

pdf('report/dataView/features_graphs/feature_scatter0.pdf')
plot(datpfeatact9[,1],datpfeatact9[,2],xlab=colnames(datpfeatact9)[1],ylab=colnames(datpfeatact9)[2],type='n',
     xlim=c(-5,(max(datpfeatact9[,1],datpfeatact9[,2])+10)))

for(i in 1:length(treatkp)){
points(datpfeatact9[datpfeatact9$Treatment %in% treatkp[i],1],
      datpfeatact9[datpfeatact9$Treatment %in% treatkp[i],2],col=i,pch=19)
  legend(25,2+0.75*i,legend=treatkp[i],col=i,bty='n',pch=19,cex=0.8)
}

dev.off()


#select those treatments with just 9 replicates to work with (treatments sent by Steffen for a start)
# keep 10 best features for now because otherwise with only 18 pairs, we'll not be estimate 

datpfeatact90<-datpfeatact[Treatment %in% treatkp,]
datpfeatact90<-datpfeatact90[,c(20,22,21,1:10)]

fit_i<-datpfeatact90[datpfeatact9$Treatment %in% treatkp[1:2],]
Y<-as.matrix(cbind(fit_i[,4:13]),ncol=10,nrow=126,byrow=T)
treat_i<-as.factor(fit_i$Treatment)

#running MANOVA
(fit <- manova(Y ~ treat_i))
summary(fit, test = "Wilks")

# https://cran.r-project.org/web/packages/HSAUR/vignettes/Ch_analysis_of_variance.pdf
summary(fit, test = "Hotelling-Lawley")
summary(fit, test = "Hotelling-Lawley")$stats[1,2]

#pvalue
summary(fit, test = "Hotelling-Lawley")$stats[1,6]
# for all the treatments in the subset, pick out the Hotelling's for plotting

#ifs
results<-array(NA,dim=c(length(treatkp),length(treatkp)))
for (j in 1:length(treatkp)) {
  for (i in 1:length(treatkp)) {
    if(j<i) { #you may as well do for (j!=i) though it is just a mirror-image
      fit_ij<-datpfeatact90[datpfeatact9$Treatment %in% treatkp[c(j,i)],]
      Y<-as.matrix(cbind(fit_ij[,4:13]),ncol=10,nrow=126,byrow=T)
      treat_i<-as.factor(fit_ij$Treatment)
      #running MANOVA
      (fit_ij <- manova(Y ~ treat_i))
      results[j,i]<-summary(fit_ij, test = "Hotelling-Lawley")$stats[1,2]
      print(summary(fit_ij, test = "Hotelling-Lawley")$stats[1,2])
    }
    else {print(paste('i=',i,'j=',j,'already done/no need to do'))}
  }
}

#distribution of pvalues
results_pvalue<-array(NA,dim=c(length(treatkp),length(treatkp)))
for (j in 1:length(treatkp)) {
  for (i in 1:length(treatkp)) {
    if(j<i) { #you may as well do for (j!=i) though it is just a mirror-image
      fit_ijpval<-datpfeatact90[datpfeatact9$Treatment %in% treatkp[c(j,i)],]
      Y<-as.matrix(cbind(fit_ijpval[,4:13]),ncol=10,nrow=126,byrow=T)
      treat_i<-as.factor(fit_ijpval$Treatment)
      #running MANOVA
      (fit_ijpval <- manova(Y ~ treat_i))
      results_pvalue[j,i]<-summary(fit_ijpval, test = "Hotelling-Lawley")$stats[1,6]
      print(summary(fit_ijpval, test = "Hotelling-Lawley")$stats[1,2])
    }
    else {print(paste('i=',i,'j=',j,'already done/no need to do'))}
  }
}

# Image of the calculated T^2 stats
# How to plot a table of numbers such that the values are represented by color?
pdf('report/dataView/manova plots/imageT2.pdf')
image(results,main='Untransformed data')
dev.off()

# into one line for plotting
results<-melt(results)
results<-results[complete.cases(results[['value']]),]
plot(results$value)

hist(results$value)


myhist <- hist(results$value)
multiplier <- myhist$counts / myhist$density
mydensity <- density(results$value)
mydensity$y <- mydensity$y * multiplier[1]

pdf('report/dataView/manova plots/hotellings_14_selectedTreatments.pdf')
plot(myhist,xlab=expression(paste('T'^2,' ','value')),ylab='counts',xlim=c(0,2000),main=paste('calculated Hotelling',"'",
                          's statistic- Untransformed data'))
lines(mydensity,lwd=2)
dev.off()
plot(mydensity,lwd=2)

#distribution of pvalues
# into one line for plotting
results_pvalue<-melt(results_pvalue)
results_pvalue<-results_pvalue[complete.cases(results_pvalue[['value']]),]

myhist_pvalue <- hist(results_pvalue$value)
multiplier_pvalue <- myhist_pvalue$counts / myhist_pvalue$density
mydensity_pvalue <- density(results_pvalue$value)
mydensity_pvalue$y <- mydensity_pvalue$y * multiplier_pvalue[1]

plot(myhist_pvalue,xlab=expression(paste('P-',' ','value')),ylab='counts',xlim=c(0,0.25))
#lines(mydensity_pvalue,lwd=2)

## yoteyote yenye ina 9 replicates??
# add the plot with lambda=10

datpfeatact10<-'p1601xxPPS_means_preprocessedfeat_actlog10.Rdata'
datpfeatact10<-file.path(datawd,datpfeatact10)
load(datpfeatact10)

#remove the endosome thing in the treatment name
datpfeatact10$Treatment=   gsub("in Colon-CytoSkeleton-Endosomes","",datpfeatact10$Treatment)

# as sent
# treatset<-c('cmpd291 @ 3 uM','cmpd788 @ 3 uM','cmpd653 @ 0.334 uM',
#             'cmpd626 @ 3 uM','cmpd104 @ 3 uM','cmpd262 @ 3 uM','cmpd298 @ 3 uM',
#             'cmpd595 @ 1 uM','cmpd464 @ 9 uM','cmpd275 @ 1 uM','cmpd306 @ 1 uM',
#             'cmpd306 @ 3 uM','cmpd486 @ 1 uM','cmpd486 @ 3 uM')

#keep relevant cmpds and concentration levels
datpfeatact10<-datpfeatact10[CMPD_ID %in% cmpd,]

#treats to keep all are 9 replicates
datpfeatact10<-datpfeatact10[Treatment %in% treatkp,]

#select those treatments with just 9 replicates to work with (treatments sent by Steffen for a start)
# keep 10 best features for now because otherwise with only 18 pairs, we'll not be estimate 

datpfeatact10<-datpfeatact10[,c(20,22,21,1:10)]


# for all the treatments in the subset, pick out the Hotelling's for plotting

#ifs
results10<-array(NA,dim=c(length(treatkp),length(treatkp)))
for (j in 1:length(treatkp)) {
  for (i in 1:length(treatkp)) {
    if(j<i) {
      fit_ij<-datpfeatact10[datpfeatact10$Treatment %in% treatkp[c(j,i)],]
      Y<-as.matrix(cbind(fit_ij[,4:13]),ncol=10,nrow=126,byrow=T)
      treat_i<-as.factor(fit_ij$Treatment)
      #running MANOVA
      (fit_ij <- manova(Y ~ treat_i))
      results10[j,i]<-summary(fit_ij, test = "Hotelling-Lawley")$stats[1,2]
      print(summary(fit_ij, test = "Hotelling-Lawley")$stats[1,2])
    }
    else {print(paste('i=',i,'j=',j,'already done/no need to do'))}
  }
}

#pvalues-distribution
results10_pvalue<-array(NA,dim=c(length(treatkp),length(treatkp)))
for (j in 1:length(treatkp)) {
  for (i in 1:length(treatkp)) {
    if(j<i) {
      fit_ij_pvalue<-datpfeatact10[datpfeatact10$Treatment %in% treatkp[c(j,i)],]
      Y<-as.matrix(cbind(fit_ij_pvalue[,4:13]),ncol=10,nrow=126,byrow=T)
      treat_i<-as.factor(fit_ij_pvalue$Treatment)
      #running MANOVA
      (fit_ij_pvalue <- manova(Y ~ treat_i))
      results10_pvalue[j,i]<-summary(fit_ij_pvalue, test = "Hotelling-Lawley")$stats[1,6]
      print(summary(fit_ij_pvalue, test = "Hotelling-Lawley")$stats[1,2])
    }
    else {print(paste('i=',i,'j=',j,'already done/no need to do'))}
  }
}

#distribution of pvalues
# into one line for plotting
results10_pvalue<-melt(results10_pvalue)
results10_pvalue<-results10_pvalue[complete.cases(results10_pvalue[['value']]),]

myhist_10pvalue <- hist(results10_pvalue$value)
multiplier_10pvalue <- myhist_10pvalue$counts / myhist_10pvalue$density
mydensity_10pvalue <- density(results10_pvalue$value)
mydensity_10pvalue$y <- mydensity_10pvalue$y * multiplier_10pvalue[1]

plot(myhist_10pvalue,xlab=expression(paste('P-',' ','value')),ylab='counts',xlim=c(0,0.25))
#lines(mydensity_pvalue,lwd=2)


# ``````````````` plots`
# Image of the calculated T^2 stats
pdf('report/dataView/manova plots/imageT2_lambda10.pdf')
image(results10,main=expression(paste(lambda,' ','=10')))
dev.off()

# into one line for plotting
results10<-melt(results10)
results10<-results10[complete.cases(results10[['value']]),]

myhist10 <- hist(results10$value)
multiplier10 <- myhist10$counts / myhist10$density
mydensity10 <- density(results10$value)
mydensity10$y <- mydensity10$y * multiplier10[1]

pdf('report/dataView/manova plots/hotellings_14_selectedTreatments_lambda10.pdf')
plot(myhist10,xlab=expression(paste('T'^2,' ','value')),ylab='counts',ylim=c(0,120),xlim=c(0,6000),
     main=expression(paste('calculated Hotelling',"'",'s statistic- transformed data',' ',
                           lambda,'=10')))
lines(mydensity10,lwd=2)
dev.off()

#lambda=25
datpfeatact25<-'p1601xxPPS_means_preprocessedfeat_actlog25.Rdata'
datpfeatact25<-file.path(datawd,datpfeatact25)
load(datpfeatact25)
datpfeatact25<-datpfeatact10
rm(datpfeatact10)
#remove the endosome thing in the treatment name
datpfeatact25$Treatment=   gsub("in Colon-CytoSkeleton-Endosomes","",datpfeatact25$Treatment)

# as sent
# treatset<-c('cmpd291 @ 3 uM','cmpd788 @ 3 uM','cmpd653 @ 0.334 uM',
#             'cmpd626 @ 3 uM','cmpd104 @ 3 uM','cmpd262 @ 3 uM','cmpd298 @ 3 uM',
#             'cmpd595 @ 1 uM','cmpd464 @ 9 uM','cmpd275 @ 1 uM','cmpd306 @ 1 uM',
#             'cmpd306 @ 3 uM','cmpd486 @ 1 uM','cmpd486 @ 3 uM')

#keep relevant cmpds and concentration levels
datpfeatact25<-datpfeatact25[CMPD_ID %in% cmpd,]

#treats to keep all are 9 replicates
datpfeatact25<-datpfeatact25[Treatment %in% treatkp,]

#select those treatments with just 9 replicates to work with (treatments sent by Steffen for a start)
# keep 10 best features for now because otherwise with only 18 pairs, we'll not be estimate 

datpfeatact25<-datpfeatact25[,c(20,22,21,1:10)]


# for all the treatments in the subset, pick out the Hotelling's for plotting

#ifs
results25<-array(NA,dim=c(length(treatkp),length(treatkp)))
for (j in 1:length(treatkp)) {
  for (i in 1:length(treatkp)) {
    if(j<i) {
      fit_ij<-datpfeatact25[datpfeatact25$Treatment %in% treatkp[c(j,i)],]
      Y<-as.matrix(cbind(fit_ij[,4:13]),ncol=25,nrow=126,byrow=T)
      treat_i<-as.factor(fit_ij$Treatment)
      #running MANOVA
      (fit_ij <- manova(Y ~ treat_i))
      results25[j,i]<-summary(fit_ij, test = "Hotelling-Lawley")$stats[1,2]
      print(summary(fit_ij, test = "Hotelling-Lawley")$stats[1,2])
    }
    else {print(paste('i=',i,'j=',j,'already done/no need to do'))}
  }
}

#pvalues-distribution
results25_pvalue<-array(NA,dim=c(length(treatkp),length(treatkp)))
for (j in 1:length(treatkp)) {
  for (i in 1:length(treatkp)) {
    if(j<i) {
      fit_ij_pvalue<-datpfeatact25[datpfeatact25$Treatment %in% treatkp[c(j,i)],]
      Y<-as.matrix(cbind(fit_ij_pvalue[,4:13]),ncol=25,nrow=126,byrow=T)
      treat_i<-as.factor(fit_ij_pvalue$Treatment)
      #running MANOVA
      (fit_ij_pvalue <- manova(Y ~ treat_i))
      results25_pvalue[j,i]<-summary(fit_ij_pvalue, test = "Hotelling-Lawley")$stats[1,6]
      print(summary(fit_ij_pvalue, test = "Hotelling-Lawley")$stats[1,2])
    }
    else {print(paste('i=',i,'j=',j,'already done/no need to do'))}
  }
}

#distribution of pvalues
# into one line for plotting
results25_pvalue<-melt(results25_pvalue)
results25_pvalue<-results25_pvalue[complete.cases(results25_pvalue[['value']]),]

myhist_25pvalue <- hist(results25_pvalue$value)
multiplier_25pvalue <- myhist_25pvalue$counts / myhist_25pvalue$density
mydensity_25pvalue <- density(results25_pvalue$value)
mydensity_25pvalue$y <- mydensity_25pvalue$y * multiplier_25pvalue[1]

plot(myhist_25pvalue,xlab=expression(paste('P-',' ','value')),ylab='counts',main=expression(paste(lambda,' ','=25')),xlim=c(0,0.25))
#lines(mydensity_pvalue,lwd=2)

pdf('report/dataView/manova plots/hotellings_14_pvalue.pdf')
par(mfrow=c(2,2))
plot(myhist_pvalue,xlab='',ylab='counts',main=expression(paste('untransformed')),xlim=c(0,0.25))
plot(myhist_10pvalue,xlab=expression(paste('P-',' ','value')),ylab='',main=expression(paste(lambda,' ','=10')),xlim=c(0,0.4))
plot(myhist_25pvalue,xlab=expression(paste('P-',' ','value')),ylab='counts',main=expression(paste(lambda,' ','=25')),xlim=c(0,0.4))
dev.off()

#scatters of P-values

pdf('report/dataView/manova plots/hotellings_14_pvalues_scatter.pdf')
plot(results_pvalue$value,results10_pvalue$value,pch=19,col='black',xlab='untransformed p-values',ylab='transformed p-values')
abline(a=0,b=1,lwd=2,col='black')
points(results_pvalue$value,results25_pvalue$value,pch=19,col='red')
legend(0.05,0.3125,legend = c(expression(paste(lambda,' ','=10')),expression(paste(lambda,' ','=25')),'equal p-value line'), bty="o",
       lty=c(0,0,1), pch=c(19,19,NA),col=c('black','red','black'),lwd=c(0,0,2))
       
dev.off()       
  
#logs are easier to see?
pdf('report/dataView/manova plots/hotellings_14_logpvalues_scatter.pdf')
plot(log(results_pvalue$value),log(results10_pvalue$value),pch=19,col='black',xlab='untransformed p-values',ylab='transformed p-values',ylim=c(-30,5),
     xlim=c(-22,0))
abline(a=0,b=1,lwd=2,col='black')
points(log(results_pvalue$value),log(results25_pvalue$value),pch=19,col='red')
legend('bottomright',legend = c(expression(paste(lambda,' ','=10')),expression(paste(lambda,' ','=25')),'equal p-value line'), bty="o",
       lty=c(0,0,1), pch=c(19,19,NA),col=c('black','red','black'),lwd=c(0,0,2))
dev.off()
       

# ``````````````` plots`
# Image of the calculated T^2 stats
pdf('report/dataView/manova plots/imageT2_lambda25.pdf')
image(results25,main=expression(paste(lambda,' ','=25')))
dev.off()

# into one line for plotting
results25<-melt(results25)
results25<-results25[complete.cases(results25[['value']]),]

myhist25 <- hist(results25$value)
multiplier25 <- myhist25$counts / myhist25$density
mydensity25 <- density(results25$value)
mydensity25$y <- mydensity25$y * multiplier25[1]

pdf('report/dataView/manova plots/hotellings_14_selectedTreatments_lambda25.pdf')
plot(myhist25,xlab=expression(paste('T'^2,' ','value')),ylab='counts',ylim=c(0,120),xlim=c(0,3500),
     main=expression(paste('calculated Hotelling',"'",'s statistic- transformed data',' ',
                           lambda,'=25')))
lines(mydensity25,lwd=2)
dev.off()

#superimpossing the densities
pdf('report/dataView/manova plots/hotellings_14_selectedTreatments_superimposed_lines.pdf')
plot(mydensity10,lwd=2,col='black',xlab=expression(paste('T'^2,' ','value')),main='')
lines(mydensity25,lwd=2,col='red',lty=1)
lines(mydensity,lwd=2,col='black',lty=2)
legend('topright',legend=c(expression(paste(lambda,' ','=10')),expression(paste(lambda,' ','=25')),'untransformed'),lty=c(1,1,2),
       col=c('black','red','black'),lwd=2)
dev.off()
##### ````````````````````````````````````````

# the distributions in terms of mean and sd
c(mean=c(mean(results$value),mean(results10$value),mean(results25$value)),sd=c(sd(results$value),sd(results10$value),sd(results25$value)))



# wait, u'll do it later

# LDA- we still have covariates that are highly correlated, making the LDA unstable

lda<- lda(formula = Treatment ~ ., 
        data = datpfeatact9[,c(1:4,'Treatment')],prior = c(1,1)/2)
cor(datpfeatact9[,-c('Treatment')])
plot(datpfeatact9[,-c('Treatment')])
#probs
lda$prior
lda$counts
lda$means
lda$scaling #linear combination coefficients (scaling) for each 
#linear discriminant 
lda$svd #singular values (svd) that gives the ratio of the between- and 
#within-group standard deviations on the linear discriminant variables



# amount of the between-group variance that is explained by each linear 
# discriminant.
prop = lda$svd^2/sum(lda$svd^2)
prop #1st LD explains 99% of the btwn grp variability

# 
# library(corrplot)
# library(caret)
# #corrplot: the library to compute correlation matrix.
# 
# #
# corMatMy <- cor(datpfeatact9[,-c('Treatment')])
# 
# #compute the correlation matrix
# 
# corrplot(corMatMy, order = "hclust")
# #visualize the matrix, clustering features by correlation index.
# 
# # inspect and set the correlation matrix at 0.5
# highlyCor <- findCorrelation(corMatMy, 0.5)
# 
# #Apply correlation filter at 0.5,
# #then we remove all the variable correlated with more 0.5
# datMyFiltered.scale <- datpfeatact9[,-highlyCor]
# corMatMy <- cor(as.numeric(datMyFiltered.scale))
# corrplot(corMatMy, order = "hclust")

#
library(Hotelling)

### in general 

z1<-matrix(c(214.969,214.823,129.943,130.3,129.72,130.193,8.305,10.53,10.168,11.133,141.517,139.45),nrow=6,byrow = T); z1

S1=matrix(c(0.150, 0.058,0.057,0.057,0.014,0.005,0.058,0.133,0.086,0.057,0.049,-0.043,0.057,0.086,0.126,0.058,0.031,
            -0.024,0.057,0.057,0.058,0.413,-0.263,-0.000,0.014,0.049,0.031,-0.263,0.421,-0.075,0.005,-0.043,-0.024,-0.000,
            -0.075,0.200),nrow=6,byrow = T)
zdiff<-matrix(z1[,1]-z1[,2],nrow = 6)

S2<-matrix(c(0.124,0.032,0.024,-0.101,0.019,0.012,0.032,0.065,0.047,-0.024,-0.012,-0.005,0.024,0.047,0.089,-0.019,0.000,0.034,-0.101,-0.024,
             -0.019,1.281,-0.490,0.238,0.019,-0.012,0.000,-0.490,0.404,-0.022,0.012,-0.005,0.034,0.238,-0.022,0.311),nrow=6,byrow=T)
Sp<-(99*S1+99*S2)/198

t(zdiff)%*% solve(Sp*(2/100))%*%zdiff

#example
######################### manual calculation
i=1; j=2
fit_ij0<-datpfeatact90[datpfeatact90$Treatment %in% treatkp[c(j,i)],]

#mean vector difference
X1_2<-matrix(apply(fit_ij0[fit_ij0$Treatment==names(table(fit_ij0$Treatment))[1],4:13],MARGIN = 2,FUN = mean),nrow=10,byrow = T)- 
  matrix(apply(fit_ij0[fit_ij0$Treatment==names(table(fit_ij0$Treatment))[2],4:13],MARGIN = 2,FUN = mean),nrow=10,byrow = T)


matrix(var(fit_ij0[fit_ij0$Treatment==names(table(fit_ij0$Treatment))[1],4:13]),ncol=10,byrow = T)
matrix(var(fit_ij0[fit_ij0$Treatment==names(table(fit_ij0$Treatment))[2],4:13]),ncol=10,byrow = T)

#pooling the variance
sp<-((as.numeric(table(fit_ij0$Treatment)[1])-1)*matrix(var(fit_ij0[
  fit_ij0$Treatment==names(table(fit_ij0$Treatment))[1],4:13]),ncol=10,byrow = T)+
    (as.numeric(table(fit_ij0$Treatment)[2])-1)*matrix(var(fit_ij0[
      fit_ij0$Treatment==names(table(fit_ij0$Treatment))[2],4:13]),ncol=10,byrow = T))/
  (as.numeric(table(fit_ij0$Treatment)[1]+as.numeric(table(fit_ij0$Treatment)[2]-2)))

#calculated hoteling
t(X1_2)%*% solve(sp*(1/as.numeric(table(fit_ij0$Treatment)[2])+ 1/as.numeric(table(fit_ij0$Treatment)[1])))%*%X1_2

# hotelling type
split.data = split(fit_ij0[,-(1:3)],fit_i$Treatment) 
x = split.data[[1]] 
y = split.data[[2]]
hotelling.stat(x, y) 
hotelling.stat(x, y, TRUE)

k=hotelling.stat(x, y)$statistic; k

###manova type
fij0<-manova(cbind(fit_ij0[,c(4:13)])~ as.factor(fit_ij0[,1]))


#