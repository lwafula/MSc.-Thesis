
# AGGREGATED DATASET2
# LWafula 31072017
# hoteling's statistic method on glog transformed data: only on active treatments
# pick the 1st 10 features [for a dataset of 18 observations, 2 df lost becoz of means vecctors for each treatment in the pair, residual rank=16
# pick number of features sufficiently lower than 16 there4 to increase power of the MANOVA-based statistic
# this results are to be used to compare with the untransformed plot for the Hotelling's results and see if there's a better separation
# alpha estimates are available from the LWafula -aggregated GLOG parameters estimation script in the MThesis/scripts folder
# the y-alpha dataset is also built in the LWafula -aggregated  GLOG parameters estimation script

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

datawd<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data'
datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy'
datawdmatrixpvals<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/results-matricesH2_pvalues etc'

# for the untransformed data is already done and saved

# add results for lambdas [experimented using lambda=2.5] then loop through all available lambdas and make a plot
file.names <- dir(datawdtranslambda, pattern ="aggregated_p1602xxPPS_means_preprocessedfeat_act_")
for (i in 1:length(file.names)){
  fileb1<-file.path(datawdtranslambda,file.names[i])
  load(fileb1)
  datpopt.featurelambda<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(glogdatp, "features_mRMR_optN")[1:10])]
  rm(list=c('glogdatp','fileb1'))
  
  
  # distribution of the Hoteling's T2 statistic and its p-value
  treatkp<-names(table(datpopt.featurelambda$Treatment))
  #lambda value
  
  lambdaval=   gsub("aggregated_p1602xxPPS_means_preprocessedfeat_act_","",file.names[i]); lambdaval<-gsub(".RData","",lambdaval); lambdaval
  resultsname<-paste0('results',lambdaval)
  results<-array(NA,dim=c(length(names(table(datpopt.featurelambda$Treatment))),length(names(table(datpopt.featurelambda$Treatment)))))
  colnames(results)<-paste(names(table(datpopt.featurelambda$Treatment)),sep = '')
  rownames(results)<-paste(names(table(datpopt.featurelambda$Treatment)),sep = '')
  
  #distribution of pvalues
  results_pvalue<-array(NA,dim=c(length(names(table(datpopt.featurelambda$Treatment))),length(names(table(datpopt.featurelambda$Treatment)))))
  colnames(results_pvalue)<-paste(names(table(datpopt.featurelambda$Treatment)),sep = '')
  rownames(results_pvalue)<-paste(names(table(datpopt.featurelambda$Treatment)),sep = '')
  
  for (j in 1:length(treatkp)) {
    for (i in 1:length(treatkp)) {
      if(j<i) { #you may as well do for (j!=i) though it is just a mirror-image
        fit_ij<-datpopt.featurelambda[datpopt.featurelambda$Treatment %in% treatkp[c(j,i)],]
        treatkp<-names(table(datpopt.featurelambda$Treatment))
        
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
  resultsname<-paste0('aggregated_results',lambdaval,'.csv')
  resultsname<-file.path(datawdmatrixpvals,resultsname)
  
  # into one line for plotting
  results<-melt(results)
  results<-results[complete.cases(results[['value']]),]
  write.csv(results,resultsname)
  
  resultsnamepvalues<-paste0('aggregated_results_pvalue',lambdaval,'.csv')
  resultsnamepvalues<-file.path(datawdmatrixpvals,resultsnamepvalues)
  results_pvalue<-melt(results_pvalue)
  results_pvalue<-results_pvalue[complete.cases(results_pvalue[['value']]),]
  write.csv(results_pvalue,resultsnamepvalues)
  
  # histograms
  myhist <- hist(results$value)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(results$value)
  mydensity$y <- mydensity$y * multiplier[1]
  
  h2<-paste0('report/dataView/manova plots/D2_aggregated_hotellings_allactiveTRTS_lambda=',lambdaval,'.pdf')
  
  pdf(h2)
  plot(myhist,xlab=expression(paste('T'^2,' ','value')),ylab='counts',xlim=c(0,summary(results$value)[6]),
       main=paste('calculated Hotelling',"'",'s statistic- transformed data'),ylim=c(0,1000000))
  lines(mydensity,lwd=2)
  dev.off()
  
}

# 31-07-2017: ALL THE ABOVE ARE ALREADY SAVED

# combining the untransformed, single-cell level transformation+well-level transformed dataset for T2
# untransformed
lambv<-c()
lambdas=c(0.1,0.5,8.5,10.5)
for (i in 1:length(lambdas)){
  h2alltrst10features<-read.table('report/results-matricesH2_pvalues etc/July2017results_H2_allTrts10fts_untransformed.csv',header = T,sep=',')
  h2alltrst10features<-h2alltrst10features[,-1]
  
  # into one line for plotting
  # h2alltrst10features<-melt(h2alltrst10features) # already done
  h2alltrst10features<-h2alltrst10features[complete.cases(h2alltrst10features[['value']]),]
  myhist <- hist(h2alltrst10features$value)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(h2alltrst10features$value)
  mydensity$y <- mydensity$y * multiplier[1]
  
  # #add here plots for other transformations for the report and H2 comparisons btwn untransformed and transformed
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'0'
  
  #single-cell level
  h2temp<-file.path(datawdmatrixpvals,paste0('results',lambdas[i],'.csv'))
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
  mydensitytemp$gr<-'1'
  mydensity<-data.frame(rbind(mydensity,mydensitytemp))
  
  # well-level
  #add transform distribution lines for Hotelling's: both single-cell+aggregated level for that transformation
  h2temp<-file.path(datawdmatrixpvals,paste0('aggregated_results',lambdas[i],'.csv'))
  lambv[i]<-lambdas[i];lambv[i]
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  
  #add the lines
  myhisttemp <- hist(h2temp$value,plot=FALSE)
  multiplier <- myhisttemp$counts / myhisttemp$density
  mydensitytemp <- density(h2temp$value)
  mydensitytemp$y <- mydensitytemp$y * multiplier[1]
  mydensitytemp<-data.frame(cbind(y=mydensitytemp$y,x=mydensitytemp$x))
  mydensitytemp$gr<-'2'
  mydensity<-data.frame(rbind(mydensity,mydensitytemp))
  
  #plot for this lambda
  mydensity$gr<-as.numeric(mydensity$gr)
  gr<-sort(unique(mydensity$gr))
  
  #plots
  hi<-paste0('report/dataView/manova plots/individual lambdas/D2_aggregated_hotellings_allactiveTRTS_lambda=',lambv[i],'.pdf')
  pdf(hi)
  
  plot(mydensity[mydensity$gr==gr[1],]$y~mydensity[mydensity$gr==gr[1],]$x,type='l',lwd=1,lty=2,ylim=c(0,max(mydensity$y)),
       col='black',xlab=expression(paste('T'^2,' ','value')),ylab='Density',
       main=expression(paste('Distribution of',' ','T'^2,' ','statistic with changing',' ',lambda)))
  lines(mydensity[mydensity$gr==1,]$y~mydensity[mydensity$gr==1,]$x,col='black',lwd=1,lty=1) #single cell level
  lines(mydensity[mydensity$gr==2,]$y~mydensity[mydensity$gr==2,]$x,col='red',lwd=1,lty=1) #welllevel
  legend('topright',legend = 
           c(expression(paste(lambda,' ','= 0 [Untransformed]')),eval(bquote(expression(paste(lambda ~'='~.(lambv[i])~'[cell level]')))),
             eval(bquote(expression(paste(lambda ~'='~.(lambv[i])~'[well level]'))))),col=c('black','black','red'),bty = 'n',cex=0.85,lwd=c(1,1,1),lty=c(2,1,1))
  dev.off()
}


############ logs for the above ############

# combining the untransformed, single-cell level transformation+well-level transformed dataset for T2
# untransformed
lambv<-c()
lambdas=c(0.1,0.5,8.5,10.5)

for(i in 1:length(lambdas)){
  myhistlog <- hist(log(h2alltrst10features$value), density = NULL)
  multiplierlog <- myhistlog$counts
  mydensitylog <- density(log(h2alltrst10features$value))
  mydensitylog$y <- mydensitylog$y * multiplierlog[1]
  
  mydensitylog<-data.frame(cbind(y=mydensitylog$y,x=mydensitylog$x))
  mydensitylog$gr<-'0'
  
  # transformed at single-cell level
  h2temp<-file.path(datawdmatrixpvals,paste0('results',lambdas[i],'.csv'))
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
  mydensitytemp$gr<-'1'
  mydensitylog<-data.frame(rbind(mydensitylog,mydensitytemp))
  
  #well-level
  h2temp<-file.path(datawdmatrixpvals,paste0('aggregated_results',lambdas[i],'.csv'))
  lambv[i]<-lambdas[i];lambv[i]
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  
  #add the lines
  myhisttemp <- hist(log(h2temp$value),plot=FALSE)
  multiplier <- myhisttemp$counts
  mydensitytemp <- density(log(h2temp$value))
  mydensitytemp$y <- mydensitytemp$y * multiplier[1]
  
  mydensitytemp<-data.frame(cbind(y=mydensitytemp$y,x=mydensitytemp$x))
  mydensitytemp$gr<-'2'
  mydensitylog<-data.frame(rbind(mydensitylog,mydensitytemp))
  
  #plot for this lambda
  mydensitylog$gr<-as.numeric(mydensitylog$gr)
  gr<-sort(unique(mydensity$gr))
  
  mydensitylog$gr<-as.numeric(mydensitylog$gr)
  gr<-sort(unique(mydensitylog$gr))
  
  # plot
  # #add here plots for other transformations for the report and H2 comparisons btwn untransformed and transformed
  hi<-paste0('report/dataView/manova plots/individual lambdas/logs/D2_aggregated_LOG_hotellings_allactiveTRTS_lambda=',lambdas[i],'.pdf')
  pdf(hi)
  
  plot(mydensitylog[mydensitylog$gr==gr[1],]$y~mydensitylog[mydensitylog$gr==gr[1],]$x,type='l',lwd=1,lty=1,ylim=c(0,max(mydensitylog$y)),xlim=c(0,14),
       col='black',xlab=expression(paste('log(',' ','T'^2,' ','value',')')),ylab='Density',
       main=expression(paste('Distribution of',' ','log(T'^2,')',' ','statistic with changing',' ',lambda)))
  lines(mydensitylog[mydensitylog$gr==1,]$y~mydensitylog[mydensitylog$gr==1,]$x,col='blue',lwd=1,lty=1) #single cell level
  lines(mydensitylog[mydensitylog$gr==2,]$y~mydensitylog[mydensitylog$gr==2,]$x,col='red',lwd=1,lty=1) #well level
  abline(v=mean(mydensitylog[mydensitylog$gr==0,]$x),lwd=1,lty=2, col='black') #mean for untransformed
  abline(v=mean(mydensitylog[mydensitylog$gr==1,]$x),lwd=1,lty=2, col='blue') #mean for single cell level
  abline(v=mean(mydensitylog[mydensitylog$gr==2,]$x),lwd=1,lty=2, col='red') #mean for well level
  
  legend('topleft',legend = 
           c(expression(paste(lambda,' ','= 0 [Untransformed]')),eval(bquote(expression(paste(lambda ~'='~.(lambv[i])~'[cell level]')))),
             eval(bquote(expression(paste(lambda ~'='~.(lambv[i])~'[well level]')))),
             expression(paste('mean untransformed')),expression(paste('mean log(','T'^2,')',' ','transformed-single cell')),
             expression(paste('mean log(','T'^2,')',' ','transformed-well level'))),col=c('black','blue','red','black','blue','red'),
         bty = 'n',cex=0.85,lwd=c(1,1,1,1,1,1), lty=c(1,1,1,2,2,2))
  
  dev.off()
}

# P-VALUES

lambv<-c()
lambdas=c(0.1,0.5,8.5,10.5)
for(i in 1:length(lambdas)){
  results_H2pvalues_allTrts10fts_untransformed<-'July2017results_H2pvalues_allTrts10fts_untransformed.csv'
  results_H2pvalues_allTrts10fts_untransformed<-file.path(datawdmatrixpvals,results_H2pvalues_allTrts10fts_untransformed)
  results_H2pvalues_allTrts10fts_untransformed<-read.table(results_H2pvalues_allTrts10fts_untransformed,header=T,sep=',')
  results_H2pvalues_allTrts10fts_untransformed<-results_H2pvalues_allTrts10fts_untransformed[,-1]
  results_H2pvalues_allTrts10fts_untransformed<-
    results_H2pvalues_allTrts10fts_untransformed[complete.cases(results_H2pvalues_allTrts10fts_untransformed$value),]
  
  results_H2pvalues_allTrts10fts_untransformed$grp<-'untransformed'
  results_H2pvalues_allTrts10fts_untransformed<-results_H2pvalues_allTrts10fts_untransformed[,-c(1,2)]
  
  #single cell level pvalues
  h2temp<-file.path(datawdmatrixpvals,paste0('results_pvalue',lambdas[i],'.csv'))
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  h2temp<-h2temp[complete.cases(h2temp$value),]
  lambv[i]<-lambdas[i];lambv[i]
  h2temp$grp<-'cell level'
  h2temp<-h2temp[,-c(1,2)]
  
  #merge for plotting and differentiating
  merge<-rbind(results_H2pvalues_allTrts10fts_untransformed,h2temp)
  
  # well level p-values
  h2temp<-file.path(datawdmatrixpvals,paste0('aggregated_results_pvalue',lambdas[i],'.csv'))
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  
  lambv[i]<-lambdas[i];lambv[i]
  h2temp$grp<-'well level'
  h2temp<-h2temp[,-c(1,2)]
  
  #merge for plotting and differentiating
  merge<-rbind(merge,h2temp)
  
  #scatter plots for each lambda against the untransformed
  hi<-paste0('report/dataView/manova plots/individual lambdas/D2_aggregated_pvalue_hotellings_allactiveTRTS_lambda=',lambv[i],'.pdf')
  pdf(hi)
  par(xpd = T, mar = par()$mar + c(0,0,0,7))
  plot(merge[['value']],type='n',xlab='index',ylab='p-value', xlim=c(0,max(length(merge[merge[['grp']]=='untransformed',][['value']]),
                                                                           length(merge[merge[['grp']]=='cell level',][['value']]),
                                                                           length(merge[merge[['grp']]=='well level',][['value']]))))
  abline(h=0.05,col='black',lwd=2,lty=1)
  points(merge[merge[['grp']]=='untransformed',][['value']],pch=19,col='black',cex=0.65)
  points(merge[merge[['grp']]=='cell level',][['value']],pch=17,col='blue',cex=0.65)
  points(merge[merge[['grp']]=='well level',][['value']],pch=15,col='red',cex=0.65)
  
  # logs? didn't gv anything good enough
  # what proportion above 0.05?
  # add legend with thz information about %
  unt_pval5<-round(length(merge[merge[['grp']]=='untransformed' & merge[['value']]>0.05,][['value']])*100/
                     length(merge[merge[['grp']]=='untransformed',][['value']]),3)
  t_pval5<-round(length(merge[merge[['grp']]=='cell level' & merge[['value']]>0.05,][['value']])*100/
                   length(merge[merge[['grp']]=='cell level',][['value']]),3)
  t0_pval5<-round(length(merge[merge[['grp']]=='well level' & merge[['value']]>0.05,][['value']])*100/
                    length(merge[merge[['grp']]=='well level',][['value']]),3)
  
  legend(max(length(merge[merge[['grp']]=='untransformed',][['value']]),length(merge[merge[['grp']]=='cell level',][['value']]),
             length(merge[merge[['grp']]=='well level',][['value']]))+
           0.1*(max(length(merge[merge[['grp']]=='untransformed',][['value']]),length(merge[merge[['grp']]=='cell level',][['value']]),
                    length(merge[merge[['grp']]=='well level',][['value']]))),0.8, 
         legend = c(paste('untransformed =',unt_pval5),paste('cell level =',t_pval5),
                    paste('well level =',t0_pval5),'0.05 significance level'),pch=c(19,17,15,NA), col=c('black','blue','red','black'), 
         lty=c(NA,NA,NA,1),lwd=c(NA,NA,NA,2),bty='n',cex=0.75,title=" % p-values>0.05")
  
  #restore mar to default value 
  par(mar=c(5, 4, 4, 2) + 0.1)      
  dev.off()
}