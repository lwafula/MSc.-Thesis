

# LWafula 14062017
# chocha ya Istanbul

# hoteling's statistic method on glog transformed data
# all active treatments
# pick the 1st 10 features [for a dataset of 18 observations, 2 df lost becoz of means for every feature, 
# 3 more due to variance-covariance? (5 df lost)-ask Ziv]

#this results are to be used to add to the untranformed plot for the Hotelling's results and see if there's a better sepration via the 
#distribution of the Hoteling's T-sqr statistic

# alpha estimates are available from the LWafula GLOG parameters estimation script in the MThesis/scripts folder
# the y-alpha dataset is also built in the LWafula GLOG parameters estimation script


# 2. includes results from lambdas [9,9.5, and 10.5 T2 statistics compared to untransformed -with related treatments identified]

# 
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

# for the untransformed data
fn_out0          <- 'p1601xxPPS_means_preprocessedfeat_act.Rdata'
fn_out0          <- file.path(datawd, fn_out0)
load(fn_out0)

# 707 treatments, top 10 features
# for every we need 18 obs (18 for each treatment)
load('data/p1601xxPPS_means_preprocessed.Rdata')
datpOpt.feature<-datp[datp[['Active_FractionOfRepl']]>0.5 ,c("Treatment", attr(datp, "features_mRMR_optN")[1:10])]
rm(list=c('datp','fn_out0','datpfeatact'))


# distribution of the Hoteling's T2 statistic and its p-value
treatkp<-names(table(datpOpt.feature$Treatment))
results<-array(NA,dim=c(length(names(table(datpOpt.feature$Treatment))),length(names(table(datpOpt.feature$Treatment)))))

#distribution of pvalues
results_pvalue<-array(NA,dim=c(length(names(table(datpOpt.feature$Treatment))),length(names(table(datpOpt.feature$Treatment)))))

for (j in 1:length(treatkp)) {
  for (i in 1:length(treatkp)) {
    if(j<i) { #you may as well do for (j!=i) though it is just a mirror-image
      fit_ij<-datpOpt.feature[datpOpt.feature$Treatment %in% treatkp[c(j,i)],]
      Y<-as.matrix(cbind(fit_ij[,2:11]),ncol=10,nrow=dim(datpOpt.feature)[1],byrow=T)
      treat_i<-as.factor(fit_ij$Treatment)
      #running MANOVA
      (fit_ij <- manova(Y ~ treat_i))
      results[j,i]<-summary(fit_ij, test = "Hotelling-Lawley")$stats[1,2]
      results_pvalue[j,i]<-summary(fit_ij, test = "Hotelling-Lawley")$stats[1,6]
      print(summary(fit_ij, test = "Hotelling-Lawley")$stats[1,2])
    }
    else {print(paste('i=',i,'j=',j,'already done/no need to do'))}
  }
}

#save the matrices b4 proceeding
# results_H2_allTrts10fts_untransformed<-'report/results-matricesH2_pvalues etc/July2017results_H2_allTrts10fts_untransformed.csv'
# write.csv(results,results_H2_allTrts10fts_untransformed)
# 
# r<-read.table('report/results-matricesH2_pvalues etc/results_H2_allTrts10fts_untransformed.csv',header = T,sep=',')
# r<-r[,-1]
# results_H2pvalues_allTrts10fts_untransformed<-'report/results-matricesH2_pvalues etc/July2017results_H2pvalues_allTrts10fts_untransformed.csv'
# write.csv(results_pvalue,results_H2pvalues_allTrts10fts_untransformed)

# hoteling's dist. plots
# How to plot a table of numbers such that the values are represented by color?
pdf('report/dataView/manova plots/imageT2_alltreats_untransformed.pdf')
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

pdf('report/dataView/manova plots/hotellings_allactiveTRTS_untransformed.pdf')
plot(myhist,xlab=expression(paste('T'^2,' ','value')),ylab='counts',xlim=c(0,22000),main=paste('calculated Hotelling',"'",
's statistic- Untransformed data'),ylim=c(0,1000000))
lines(mydensity,lwd=2)
dev.off()

#add here plots for other transformations for the report and H2 comparisons btwn untransformed and transformed
plot(mydensity,lwd=2,col='black',xlab=expression(paste('T'^2,' ','value')),main='')


# pvalues
# into one line for plotting
results_pvalue<-melt(results_pvalue)
results_pvalue<-results_pvalue[complete.cases(results_pvalue[['value']]),]

myhist_pvalue <- hist(results_pvalue$value)
multiplier_pvalue <- myhist_pvalue$counts / myhist_pvalue$density
mydensity_pvalue <- density(results_pvalue$value)
mydensity_pvalue$y <- mydensity_pvalue$y * multiplier_pvalue[1]

plot(myhist_pvalue,xlab=expression(paste('P-',' ','value')),ylab='counts',xlim=c(0,0.25))

# logs and add other transformations pvalues
  # plot(log(results_pvalue$value),log(results10_pvalue$value), pch=19,col='black',xlab='untransformed p-values',ylab='transformed p-values',ylim=c(-30,5),
  #      xlim=c(-22,0))
  # abline(a=0,b=1,lwd=2,col='black')
  # points(log(results_pvalue$value),log(results25_pvalue$value),pch=19,col='red')
  # legend('bottomright',legend = c(expression(paste(lambda,' ','=10')),expression(paste(lambda,' ','=25')),'equal p-value line'), bty="o",
  #        lty=c(0,0,1), pch=c(19,19,NA),col=c('black','red','black'),lwd=c(0,0,2))
  # 


# add results for lambdas [experimented using lambda=2.5] then loop through all available lambdas and make a plot
file.names <- dir(datawdtranslambda, pattern ="p1601xxPPS_means_preprocessedfeat_act_")
for (i in 1:length(file.names)){
fileb1<-file.path(datawdtranslambda,file.names[i])
load(fileb1)
datpopt.featurelambda<-glogdatp[glogdatp[['Active_FractionOfRepl']]>0.5 ,c("Treatment", attr(glogdatp, "features_mRMR_optN")[1:10])]
rm(list=c('glogdatp','fileb1'))


# distribution of the Hoteling's T2 statistic and its p-value
treatkp<-names(table(datpopt.featurelambda$Treatment))
#lambda value

lambdaval=   gsub("p1601xxPPS_means_preprocessedfeat_act_","",file.names[i]); lambdaval<-gsub(".RData","",lambdaval); lambdaval
resultsname<-paste0('results',lambdaval)
results<-array(NA,dim=c(length(names(table(datpopt.featurelambda$Treatment))),length(names(table(datpopt.featurelambda$Treatment)))))

#distribution of pvalues
results_pvalue<-array(NA,dim=c(length(names(table(datpopt.featurelambda$Treatment))),length(names(table(datpopt.featurelambda$Treatment)))))

for (j in 1:length(treatkp)) {
  for (i in 1:length(treatkp)) {
    if(j<i) { #you may as well do for (j!=i) though it is just a mirror-image
      fit_ij<-datpopt.featurelambda[datpopt.featurelambda$Treatment %in% treatkp[c(j,i)],]
      Y<-as.matrix(cbind(fit_ij[,2:11]),ncol=10,nrow=dim(datpopt.featurelambda)[1],byrow=T)
      treat_i<-as.factor(fit_ij$Treatment)
      #running MANOVA
      (fit_ij <- manova(Y ~ treat_i))
      results[j,i]<-summary(fit_ij, test = "Hotelling-Lawley")$stats[1,2]
      results_pvalue[j,i]<-summary(fit_ij, test = "Hotelling-Lawley")$stats[1,6]
      print(summary(fit_ij, test = "Hotelling-Lawley")$stats[1,2])
    }
    else {print(paste('i=',i,'j=',j,'already done/no need to do'))}
  }
}

#save the matrices b4 proceeding
resultsname<-paste0('results',lambdaval,'.csv')
resultsname<-file.path(datawdmatrixpvals,resultsname)
write.csv(results,resultsname)

resultsnamepvalues<-paste0('results_pvalue',lambdaval,'.csv')
resultsnamepvalues<-file.path(datawdmatrixpvals,resultsnamepvalues)
write.csv(results_pvalue,resultsnamepvalues)

imageplot<-paste0('report/dataView/manova plots/image lambda=',lambdaval,'.pdf')
lambdaval<-as.numeric(lambdaval)
pdf(imageplot)
image(results,main=expression(paste('transformed data')))
text(0.75, 0.1, expression(lambda),cex = .9)
text(0.79, 0.1,'=')
text(0.84, 0.1, lambdaval)
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

h2<-paste0('report/dataView/manova plots/hotellings_allactiveTRTS_lambda=',lambdaval,'.pdf')

pdf(h2)
plot(myhist,xlab=expression(paste('T'^2,' ','value')),ylab='counts',xlim=c(0,10000),
     main=paste('calculated Hotelling',"'",'s statistic- transformed data'),ylim=c(0,1000000))
lines(mydensity,lwd=2)
dev.off()


# pvalues
# into one line for plotting
results_pvalue<-melt(results_pvalue)
results_pvalue<-results_pvalue[complete.cases(results_pvalue[['value']]),]

myhist_pvalue <- hist(results_pvalue$value)
multiplier_pvalue <- myhist_pvalue$counts / myhist_pvalue$density
mydensity_pvalue <- density(results_pvalue$value)
mydensity_pvalue$y <- mydensity_pvalue$y * multiplier_pvalue[1]

plot(myhist_pvalue,xlab=expression(paste('P-',' ','value')),ylab='counts',xlim=c(0,0.25))

}


# 18-06-2017: ALL THE ABOVE ARE ALREADY SAVED
# Hotelling's distribution

h2alltrst10features<-read.table('report/results-matricesH2_pvalues etc/July2017results_H2_allTrts10fts_untransformed.csv',header = T,sep=',')
h2alltrst10features<-h2alltrst10features[,-1]

# into one line for plotting
h2alltrst10features<-melt(h2alltrst10features)
h2alltrst10features<-h2alltrst10features[complete.cases(h2alltrst10features[['value']]),]

myhist <- hist(h2alltrst10features$value)
multiplier <- myhist$counts / myhist$density
mydensity <- density(h2alltrst10features$value)
mydensity$y <- mydensity$y * multiplier[1]

# #add here plots for other transformations for the report and H2 comparisons btwn untransformed and transformed
hi<-paste0('report/dataView/manova plots/individual lambdas/hotellings_allactiveTRTS_lambda=','all lambdas','.pdf')
pdf(hi)

mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
mydensity$gr<-'0'
plot(mydensity$y~mydensity$x,type='l',lwd=2.5,lty=2,ylim=c(0,70e+05),
     col='black',xlab=expression(paste('T'^2,' ','value')),ylab='Density',
     main=expression(paste('Distribution of',' ','T'^2,' ','statistic with changing',' ',lambda)))

#add transform distribution lines for Hotelling's
file.names <- dir(datawdmatrixpvals, pattern ="results")

#keep only Hotelling's results for now
file.names<-file.names[!grepl('^(results_)',file.names, ig = TRUE)]; file.names<-file.names[!grepl('^(July2017)',file.names, ig = TRUE)]

lambv<-c()
for (i in 1:length(file.names)){
  h2temp<-file.path(datawdmatrixpvals,file.names[i])
  lambv[i]<-gsub('results','',file.names[i]); lambv[i]<-gsub('.csv','',lambv[i]);lambv[i]
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  
  # into one line for plotting
  h2temp<-melt(h2temp)
  h2temp<-h2temp[complete.cases(h2temp[['value']]),]
  
  #add the lines
  myhisttemp <- hist(h2temp$value,plot=FALSE)
  multiplier <- myhisttemp$counts / myhisttemp$density
  mydensitytemp <- density(h2temp$value)
  mydensitytemp$y <- mydensitytemp$y * multiplier[1]
  
  mydensitytemp<-data.frame(cbind(y=mydensitytemp$y,x=mydensitytemp$x))
  mydensitytemp$gr<-lambv[i]
  mydensity<-data.frame(rbind(mydensity,mydensitytemp))
}
mydensity$gr<-as.numeric(mydensity$gr)
gr<-sort(unique(mydensity$gr))
for(j in 1:length(gr)){
  lines(mydensity[mydensity$gr==gr[j],]$y~mydensity[mydensity$gr==gr[j],]$x,col=j,lwd=2)
  legend(15000,(52e+05)*0.071*gr[j],legend = gr[j],col=j,lty=1,bty = 'n',cex=0.5,lwd=2)
}
dev.off()


#ones that seem betterly different- individual plots
file.names <- dir(datawdmatrixpvals, pattern ="results")
#keep only Hotelling's results for now
file.names<-file.names[!grepl('^(results_)',file.names, ig = TRUE)]

lambv<-c()

for (i in 1:length(file.names)){
  h2temp<-file.path(datawdmatrixpvals,file.names[i])
  lambv[i]<-gsub('results','',file.names[i]); lambv[i]<-gsub('.csv','',lambv[i]); lambv[i]
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  
  # into one line for plotting
  h2temp<-melt(h2temp)
  h2temp<-h2temp[complete.cases(h2temp[['value']]),]
  
  #add the lines
  myhisttemp <- hist(h2temp$value,plot=FALSE)
  multiplier <- myhisttemp$counts / myhisttemp$density
  mydensitytemp <- density(h2temp$value)
  mydensitytemp$y <- mydensitytemp$y * multiplier[1]
  
  #plots for each lambda
  hi<-paste0('report/dataView/manova plots/individual lambdas/hotellings_allactiveTRTS_lambda=',lambv[i],'.pdf')
  pdf(hi)
  #look for a way of avoiding thz plot inside here!
  plot(mydensity,lwd=2,col='black',xlab=expression(paste('T'^2,' ','value')),
       main=expression(paste('Distribution of',' ','T'^2,' ','statistic with changing',' ',lambda)),lty=2,
       ylim=c(0,round(summary(mydensitytemp$y, mydensity$y)[6],0)))
  lines(mydensitytemp,col=i,lwd=2)
  legend('topright',legend = c('untransformed',lambv[i]),col=c('black',i),lty=c(2,1),bty = 'n',cex=1.5)
  dev.off()
}

# p-values distribution

results_H2pvalues_allTrts10fts_untransformed<-'results_H2pvalues_allTrts10fts_untransformed.csv'
results_H2pvalues_allTrts10fts_untransformed<-file.path(datawdmatrixpvals,results_H2pvalues_allTrts10fts_untransformed)
results_H2pvalues_allTrts10fts_untransformed<-read.table(results_H2pvalues_allTrts10fts_untransformed,header=T,sep=',')
results_H2pvalues_allTrts10fts_untransformed<-results_H2pvalues_allTrts10fts_untransformed[,-1]

plot(log(results_H2pvalues_allTrts10fts_untransformed[['value']]),xlab=expression(paste('P-',' ','value')),ylab='counts')

results_H2pvalues_allTrts10fts_untransformed$grp<-'untransformed'
results_H2pvalues_allTrts10fts_untransformed<-results_H2pvalues_allTrts10fts_untransformed[,-c(1,2)]


# individual pvalues
file.names <- dir(datawdmatrixpvals, pattern ="results")
#keep only Hotelling's pvalues results for now
file.names<-file.names[grepl('^(results_pvalue)',file.names, ig = TRUE)]

lambv<-c()

for (i in 1:length(file.names)){
  h2temp<-file.path(datawdmatrixpvals,file.names[i])
  lambv[i]<-gsub('results_pvalue','',file.names[i]); lambv[i]<-gsub('.csv','',lambv[i]); lambv[i]
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  
  # into one line for plotting
  h2temp<-melt(h2temp)
  h2temp<-h2temp[complete.cases(h2temp[['value']]),]
  h2temp$grp<-paste('lambda',lambv[i])
  h2temp<-h2temp[,-c(1,2)]
  
  #merge for plotting and differentiating
  merge<-rbind(results_H2pvalues_allTrts10fts_untransformed,h2temp)
  #scatter plots for each lambda against the untransformed
  hi<-paste0('report/dataView/manova plots/individual lambdas/pvalue_hotellings_allactiveTRTS_lambda=',lambv[i],'.pdf')
  pdf(hi)
  par(xpd = T, mar = par()$mar + c(0,0,0,7))
  plot(merge[['value']],type='n',xlab='index',ylab='p-value',xlim=c(0,max(length(merge[merge[['grp']]=='untransformed',][['value']]),
                                                                          length(merge[merge[['grp']]==paste('lambda',lambv[i]),][['value']]))))
  abline(h=0.05,col='blue',lwd=2,lty=1)
  points(merge[merge[['grp']]=='untransformed',][['value']],pch=19,col='black',cex=0.65)
  points(merge[merge[['grp']]==paste('lambda',lambv[i]),][['value']],pch=19,col='red',cex=0.65)
  
  # logs? didn't gv anything good enough
  # what proportion above 0.05?
  unt_pval5<-round(length(merge[merge[['grp']]=='untransformed' & merge[['value']]>0.05,][['value']])*100/
                     length(merge[merge[['grp']]=='untransformed',][['value']]),3)
  t_pval5<-round(length(merge[merge[['grp']]==paste('lambda',lambv[i]) & merge[['value']]>0.05,][['value']])*100/
                   length(merge[merge[['grp']]==paste('lambda',lambv[i]),][['value']]),3)
  
  #add legend with thz information about %
  legend(260000,0.8, legend = c(
    paste('untransformed =',unt_pval5),paste('transformed =',t_pval5),'0.05 significance level'),pch=c(19,19,NA), col=c('black','red','blue'), 
    lty=c(NA,NA,1),lwd=c(NA,NA,2),bty='n',cex=0.75,title=" % p-values>0.05")
  
  #restore mar to default value 
  par(mar=c(5, 4, 4, 2) + 0.1)      
  dev.off()
}

#just scatter plots are not very informative, so instead combine them and plot the values differentiating by colour and cut-off for 
# untransformed and lambda corrected values


