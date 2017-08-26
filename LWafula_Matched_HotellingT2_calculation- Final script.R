

### comparison for Hotelling's T2 for only treatments that were active pre- and post-transformation [using lambda= 10.5]


# LWafula 11082017
# hoteling's statistic method on glog transformed data: only on active treatments
# pick the 1st 10 features [for a dataset of 18 observations, 2 df lost becoz of means vecctors for each treatment in the pair, residual rank=16
# pick number of features sufficiently lower than 16 there4 to increase power of the MANOVA-based statistic
# this results are to be used to compare with the untranformed plot for the Hotelling's results and see if there's a better separation
# alpha estimates are available from the LWafula GLOG parameters estimation script in the MThesis/scripts folder
# the y-alpha dataset is also built in the LWafula GLOG parameters estimation script

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
datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/data/glog transformed datasets by lambda - Copy'
datawdmatrixpvals<-'E:/uhasselt/22yrsem/MThesis/report/results-matricesH2_pvalues etc'

# 707 treatments, top 10 features
load('data/p1601xxPPS_means_preprocessed.Rdata')
datpOpt.feature<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features_mRMR_optN")[1:10])]
rm(list=c('datp'))

# the transformed data
file.names <- dir(datawdtranslambda, pattern ="^p1601xxPPS_means_preprocessedfeat_act_10.5")
  fileb1<-file.path(datawdtranslambda,file.names)
  load(fileb1)
  datpopt.featurelambda<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(glogdatp, "features_mRMR_optN")[1:10])]
  rm(list=c('glogdatp','fileb1'))

# for both, just keep those treatments that were active in both untransformed and transformed data sets
intersect(names(table(datpOpt.feature$Treatment)),names(table(datpopt.featurelambda$Treatment)) ) # 683

#datpOpt.feature.temp<-datpOpt.feature
datpOpt.feature<- 
  datpOpt.feature[datpOpt.feature$Treatment %in% intersect(names(table(datpOpt.feature$Treatment)),names(table(datpopt.featurelambda$Treatment))),]
datpopt.featurelambda<- 
  datpopt.featurelambda[datpopt.featurelambda$Treatment %in% intersect(names(table(datpOpt.feature$Treatment)),names(table(datpopt.featurelambda$Treatment))),]

 
# distribution of the Hoteling's T2 statistic and its p-value
treatkp<-names(table(datpOpt.feature$Treatment))
results<-array(NA,dim=c(length(names(table(datpOpt.feature$Treatment))),length(names(table(datpOpt.feature$Treatment)))))
colnames(results)<-paste(treatkp,sep = '')
rownames(results)<-paste(treatkp,sep = '')

#distribution of pvalues
results_pvalue<-array(NA,dim=c(length(names(table(datpOpt.feature$Treatment))),length(names(table(datpOpt.feature$Treatment)))))
colnames(results_pvalue)<-paste(treatkp,sep = '')
rownames(results_pvalue)<-paste(treatkp,sep = '')

# NOTE: thz takes time (2hrs), run without assigning the names if needed, otherwise use the saved results
for (j in 1:length(treatkp)) {
  for (i in 1:length(treatkp)) {
    if(j<i) { #you may as well do for (j!=i) though it is just a mirror-image
      fit_ij<-datpOpt.feature[datpOpt.feature$Treatment %in% treatkp[c(j,i)],]
      treatkp<-names(table(datpOpt.feature$Treatment))
      
      #calculating the Hotelling and pvalue
      (fit<-hotelling.test(.~Treatment, data = fit_ij[,c(1,2:11)]))
      #running MANOVA
      results[j,i]<-fit$stats$statistic
      
      results_pvalue[j,i]<-fit$pval
      print(fit$stats$statistic)
    }
    else {print(paste('i=',i,'j=',j,'already done/no need to do'))}
  }
}

#save the matrices b4 proceeding
results<-melt(results)
results<-results[complete.cases(results[['value']]),]

results_H2_allTrts10fts_untransformed<-'report/results-matricesH2_pvalues etc/matched_Aug2017results_H2_allTrts10fts_untransformed.csv'
write.csv(results,results_H2_allTrts10fts_untransformed)

# hoteling's dist. plots
# How to plot a table of numbers such that the values are represented by color?
myhist <- hist(results$value)
multiplier <- myhist$counts / myhist$density
mydensity <- density(results$value)
mydensity$y <- mydensity$y * multiplier[1]

pdf('report/dataView/manova plots/matched_hotellings_allactiveTRTS_untransformed.pdf')
plot(myhist,xlab=expression(paste('T'^2,' ','value')),ylab='counts',xlim=c(0,350000),main=paste('calculated Hotelling',"'",'s statistic- Untransformed data'),
     ylim=c(0,1000000))
lines(mydensity,lwd=2)
dev.off()

#add here plots for other transformations for the report and H2 comparisons btwn untransformed and transformed
plot(mydensity,lwd=2,col='black',xlab=expression(paste('T'^2,' ','value')),main='',xlim=c(0,350000))

# pvalues
# into one line for plotting
results_pvalue<-melt(results_pvalue)
results_pvalue<-results_pvalue[complete.cases(results_pvalue[['value']]),]
results_H2pvalues_allTrts10fts_untransformed<-'report/results-matricesH2_pvalues etc/matched_Aug2017results_H2pvalues_allTrts10fts_untransformed.csv'
write.csv(results_pvalue,results_H2pvalues_allTrts10fts_untransformed)

myhist_pvalue <- hist(results_pvalue$value)
multiplier_pvalue <- myhist_pvalue$counts / myhist_pvalue$density
mydensity_pvalue <- density(results_pvalue$value)
mydensity_pvalue$y <- mydensity_pvalue$y * multiplier_pvalue[1]

plot(myhist_pvalue,xlab=expression(paste('P-',' ','value')),ylab='counts',xlim=c(0,0.25))

# add results for lambda= 10.5   # distribution of the Hoteling's T2 statistic and its p-value
  treatkp<-names(table(datpopt.featurelambda$Treatment))
#lambda value 
  lambdaval=   gsub("p1601xxPPS_means_preprocessedfeat_act_","",file.names[1]); lambdaval<-gsub(".RData","",lambdaval); lambdaval
  resultsname<-paste0('results',lambdaval)
  results<-array(NA,dim=c(length(names(table(datpopt.featurelambda$Treatment))),length(names(table(datpopt.featurelambda$Treatment)))))
  
  colnames(results)<-paste(treatkp,sep = '')
  rownames(results)<-paste(treatkp,sep = '')
  
  #distribution of pvalues
  results_pvalue<-array(NA,dim=c(length(names(table(datpopt.featurelambda$Treatment))),length(names(table(datpopt.featurelambda$Treatment)))))
  colnames(results_pvalue)<-paste(treatkp,sep = '')
  rownames(results_pvalue)<-paste(treatkp,sep = '')
  
  for (j in 1:length(treatkp)) {
    for (i in 1:length(treatkp)) {
      if(j<i) { #you may as well do for (j!=i) though it is just a mirror-image
        fit_ij<-datpopt.featurelambda[datpopt.featurelambda$Treatment %in% treatkp[c(j,i)],]
        
        #calculating the Hotelling and pvalue
        (fit<-hotelling.test(.~Treatment, data = fit_ij[,c(1,2:11)]))
        #running MANOVA
        results[j,i]<-fit$stats$statistic
        
        results_pvalue[j,i]<-fit$pval
        print(fit$stats$statistic)
      }
      else {print(paste('i=',i,'j=',j,'already done/no need to do'))}
    }
  }
  

  #save the matrices b4 proceeding
  resultsname<-paste0('matched_results',lambdaval,'.csv')
  resultsname<-file.path(datawdmatrixpvals,resultsname)
  results<-melt(results)
  results<-results[complete.cases(results[['value']]),]
  write.csv(results,resultsname)
  
  resultsnamepvalues<-paste0('matched_results_pvalue',lambdaval,'.csv')
  resultsnamepvalues<-file.path(datawdmatrixpvals,resultsnamepvalues)
  results_pvalue<-melt(results_pvalue)
  results_pvalue<-results_pvalue[complete.cases(results_pvalue[['value']]),]
  write.csv(results_pvalue,resultsnamepvalues)

  # into one line for plotting
  myhist <- hist(results$value)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(results$value)
  mydensity$y <- mydensity$y * multiplier[1]
  
  h2<-paste0('report/dataView/manova plots/matched_hotellings_allactiveTRTS_lambda=',lambdaval,'.pdf')
  
  pdf(h2)
  plot(myhist,xlab=expression(paste('T'^2,' ','value')),ylab='counts',xlim=c(0,summary(results$value)[6]),
       main=paste('calculated Hotelling',"'",'s statistic- transformed data'),ylim=c(0,1000000))
  lines(mydensity,lwd=2)
  dev.off()

# 11-08-2017: ALL THE ABOVE ARE ALREADY SAVED
# Hotelling's distribution

h2alltrst10features<-read.table('report/results-matricesH2_pvalues etc/matched_Aug2017results_H2_allTrts10fts_untransformed.csv',header = T,sep=',')
h2alltrst10features<-h2alltrst10features[,-1]

# into one line for plotting
h2alltrst10features<-h2alltrst10features[complete.cases(h2alltrst10features[['value']]),]

myhist <- hist(h2alltrst10features$value)
multiplier <- myhist$counts / myhist$density
mydensity <- density(h2alltrst10features$value)
mydensity$y <- mydensity$y * multiplier[1]

# #add here plots for other transformations for the report and H2 comparisons btwn untransformed and transformed
mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
mydensity$gr<-'0'

#add transform distribution lines for Hotelling's
file.names <- dir(datawdmatrixpvals, pattern ="matched_results")

#keep only Hotelling's results for now
file.names<-file.names[!grepl('^(matched_results_)',file.names, ig = TRUE)]; file.names<-file.names[!grepl('^(Aug2017)',file.names, ig = TRUE)]

lambv<-c()
  h2temp<-file.path(datawdmatrixpvals,file.names[1])
  lambv[1]<-gsub('matched_results','',file.names[1]); lambv[1]<-gsub('.csv','',lambv[1]);lambv[1]
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  
  # into one line for plotting
  h2temp<-h2temp[complete.cases(h2temp[['value']]),]
  
  #add the lines
  myhisttemp <- hist(h2temp$value,plot=FALSE)
  multiplier <- myhisttemp$counts / myhisttemp$density
  mydensitytemp <- density(h2temp$value)
  mydensitytemp$y <- mydensitytemp$y * multiplier[1]
  
  mydensitytemp<-data.frame(cbind(y=mydensitytemp$y,x=mydensitytemp$x))
  mydensitytemp$gr<-lambv[1]
  mydensity<-data.frame(rbind(mydensity,mydensitytemp))

mydensity$gr<-as.numeric(mydensity$gr)
gr<-sort(unique(mydensity$gr))

#plots
hi<-paste0('report/dataView/manova plots/individual lambdas/matched_hotellings_allactiveTRTS_lambda=','all lambdas','.pdf')
pdf(hi)

plot(mydensity[mydensity$gr==gr[1],]$y~mydensity[mydensity$gr==gr[1],]$x,type='l',lwd=1,lty=2,ylim=c(0,45e+05),
     col='black',xlab=expression(paste('T'^2,' ','value')),ylab='Density',
     main=expression(paste('Distribution of',' ','T'^2,' ','statistic with changing',' ',lambda)))

# for(j in 1:length(gr)){
    lines(mydensity[mydensity$gr==gr[2],]$y~mydensity[mydensity$gr==gr[2],]$x,col='black',lwd=1,lty=1)
    legend('topright',legend = c(expression(paste(lambda,' ','= 0 [Untransformed]')), 
                            eval(bquote(expression(lambda~'='  ~.(lambv[1]))))),col='black',lty=c(2,1),bty = 'n',cex=0.5,lwd=1) 
# }
dev.off()

############ logs for the above ############

# logs might be better

myhistlog <- hist(log(h2alltrst10features$value), density = NULL)
multiplierlog <- myhistlog$counts
mydensitylog <- density(log(h2alltrst10features$value))
mydensitylog$y <- mydensitylog$y * multiplierlog[1]

mydensitylog<-data.frame(cbind(y=mydensitylog$y,x=mydensitylog$x))
mydensitylog$gr<-'0'

#lambda results
file.names <- dir(datawdmatrixpvals, pattern ="matched_results")

#keep only Hotelling's results for now
file.names<-file.names[!grepl('^(matched_results_pvalue10.5)',file.names, ig = TRUE)]

lambv<-c()
# for (i in 1:length(file.names)){
  h2temp<-file.path(datawdmatrixpvals,file.names[1])
  lambv[1]<-gsub('matched_results','',file.names[1]); lambv[1]<-gsub('.csv','',lambv[1]);lambv[1]
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  
  # into one line for plotting
  h2temp<-h2temp[complete.cases(h2temp[['value']]),]
  
  #add the lines
  myhisttemp <- hist(log(h2temp$value),plot=FALSE)
  multiplier <- myhisttemp$counts
  mydensitytemp <- density(log(h2temp$value))
  mydensitytemp$y <- mydensitytemp$y * multiplier[1]
  
  mydensitytemp<-data.frame(cbind(y=mydensitytemp$y,x=mydensitytemp$x))
  mydensitytemp$gr<-lambv[1]
  mydensitylog<-data.frame(rbind(mydensitylog,mydensitytemp))
# }
mydensitylog$gr<-as.numeric(mydensitylog$gr)
gr<-sort(unique(mydensitylog$gr))

# plot
# #add here plots for other transformations for the report and H2 comparisons btwn untransformed and transformed
hi<-paste0('report/dataView/manova plots/individual lambdas/logs/matched_LOG_hotellings_allactiveTRTS_lambda=',lambv[1],'.pdf')
pdf(hi)

plot(mydensitylog[mydensitylog$gr==gr[1],]$y~mydensitylog[mydensitylog$gr==gr[1],]$x,type='l',lwd=2,lty=1,ylim=c(0,2),xlim=c(0,14),
     col='black',xlab=expression(paste('log(',' ','T'^2,' ','value',')')),ylab='Density',
     main=expression(paste('Distribution of',' ','log(T'^2,')',' ','for matching active treatments')))

# for(j in 1:length(gr)){
#   if(gr[j]<1){
    lines(mydensitylog[mydensitylog$gr==gr[2],]$y~mydensitylog[mydensitylog$gr==gr[2],]$x,col='blue',lwd=2,lty=1)
    abline(v=mean(mydensitylog[mydensitylog$gr==0,]$x),lwd=1,lty=2, col='black')
    abline(v=mean(mydensitylog[mydensitylog$gr==lambv[1],]$x),lwd=1,lty=2, col='blue')
    
    legend('topleft',
           legend = c('untransformed',eval(bquote(expression(lambda~'='  ~.(lambv[1])))), expression(paste('mean log(','T'^2,')',' ','untransformed')),expression(paste('mean log(','T'^2,')',' ','transformed'))),
           col=c('black','blue','black','blue'),lty=c(1,1,2,2),bty = 'n',cex=1.5)
  # }
    # }
dev.off()

# P-VALUES

results_H2pvalues_allTrts10fts_untransformed<-'matched_Aug2017results_H2pvalues_allTrts10fts_untransformed.csv'
results_H2pvalues_allTrts10fts_untransformed<-file.path(datawdmatrixpvals,results_H2pvalues_allTrts10fts_untransformed)
results_H2pvalues_allTrts10fts_untransformed<-read.table(results_H2pvalues_allTrts10fts_untransformed,header=T,sep=',')
results_H2pvalues_allTrts10fts_untransformed<-results_H2pvalues_allTrts10fts_untransformed[,-1]

results_H2pvalues_allTrts10fts_untransformed<-results_H2pvalues_allTrts10fts_untransformed[complete.cases(results_H2pvalues_allTrts10fts_untransformed$value),]
results_H2pvalues_allTrts10fts_untransformed$grp<-'untransformed'
results_H2pvalues_allTrts10fts_untransformed<-results_H2pvalues_allTrts10fts_untransformed[,-c(1,2)]


# individual pvalues
file.names <- dir(datawdmatrixpvals, pattern ="matched_results")
#keep only Hotelling's pvalues results for now
file.names<-file.names[grepl('^(matched_results_pvalue)',file.names, ig = TRUE)]
lambv<-c()

# for (i in 1:length(file.names)){
  h2temp<-file.path(datawdmatrixpvals,file.names[1])
  lambv[1]<-gsub('matched_results_pvalue','',file.names[1]); lambv[1]<-gsub('.csv','',lambv[1]); lambv[1]
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  
  # into one line for plotting
  h2temp<-h2temp[complete.cases(h2temp[['value']]),]
  h2temp$grp<-paste('lambda',lambv[1])
  h2temp<-h2temp[,-c(1,2)]
  
  #merge for plotting and differentiating
  merge<-rbind(results_H2pvalues_allTrts10fts_untransformed,h2temp)
  #scatter plots for each lambda against the untransformed
  hi<-paste0('report/dataView/manova plots/individual lambdas/matched_pvalue_hotellings_allactiveTRTS_lambda=',lambv[1],'.pdf')
  pdf(hi)
  par(xpd = T, mar = par()$mar + c(0,0,0,7))
  plot(merge[['value']],type='n',xlab='index',ylab='p-value',xlim=c(0,max(length(merge[merge[['grp']]=='untransformed',][['value']]),
                                                                          length(merge[merge[['grp']]==paste('lambda',lambv[1]),][['value']]))))
  abline(h=0.05/(683*682*0.5),col='blue',lwd=2,lty=1)
  points(merge[merge[['grp']]=='untransformed',][['value']],pch=19,col='black',cex=0.65)
  points(merge[merge[['grp']]==paste('lambda',lambv[1]),][['value']],pch=19,col='red',cex=0.65)
  
  # logs? didn't gv anything good enough
  # what proportion above 0.05?
  # unt_pval5<-round(length(merge[merge[['grp']]=='untransformed' & merge[['value']]>0.05,][['value']])*100/
  #                    length(merge[merge[['grp']]=='untransformed',][['value']]),3)
  # t_pval5<-round(length(merge[merge[['grp']]==paste('lambda',lambv[1]) & merge[['value']]>0.05,][['value']])*100/
  #                  length(merge[merge[['grp']]==paste('lambda',lambv[1]),][['value']]),3)
  # c(unt_pval5*683*682*0.5*0.01,t_pval5*683*682*0.5*0.01,unt_pval5*683*682*0.5*0.01-t_pval5*683*682*0.5*0.01) #c( 713.9907 652.2624  61.7283)
  # 
  # proportions after correcting for multiplicity (Bonferroni)  (116020.63 106967.69   9052.94) [b4=49.815% , after=45.928%]
  unt_pval5<-round(length(merge[merge[['grp']]=='untransformed' & merge[['value']]>(0.05/(683*682*0.5)),][['value']])*100/
                     length(merge[merge[['grp']]=='untransformed',][['value']]),3)
  t_pval5<-round(length(merge[merge[['grp']]==paste('lambda',lambv[1]) & merge[['value']]>(0.05/(683*682*0.5)),][['value']])*100/
                   length(merge[merge[['grp']]==paste('lambda',lambv[1]),][['value']]),3)
  c(unt_pval5*683*682*0.5*0.01,t_pval5*683*682*0.5*0.01,unt_pval5*683*682*0.5*0.01-t_pval5*683*682*0.5*0.01) #c(305.1029-260.8514,44.2515)
  
  #add legend with thz information about %
  legend(max(length(merge[merge[['grp']]=='untransformed',][['value']]),length(merge[merge[['grp']]=='lambda 10.5',][['value']]))+
         0.035*max(length(merge[merge[['grp']]=='untransformed',][['value']]),length(merge[merge[['grp']]=='lambda 10.5',][['value']])),0.8, 
         legend = c(
    paste('untransformed =',unt_pval5),paste('transformed =',t_pval5),'significance level'),pch=c(19,19,NA), col=c('black','red','blue'), 
    lty=c(NA,NA,1),lwd=c(NA,NA,2),bty='n',cex=0.75,title=" % not significant")
  
  #restore mar to default value 
  par(mar=c(5, 4, 4, 2) + 0.1)  

  dev.off()

