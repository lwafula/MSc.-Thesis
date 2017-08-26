
# DATA SET 2
# wellvs cell vs untransformd for lambda =0.1 and 0.5
### MATCHED: comparison for Hotelling's T2 for only treatments that were active pre- and post-transformation [using lambda= 10.5]

# NOTES
# - the 1st part calculates the Hotelling's T2 for each of the lambda values using untransformed data , cell and well level data sets 
# - with the same number of actively-called treatments. This is to enable comparativity of the results
# - 2nd part plots the distributions for T2 and associated p-values. Notice that the outputs are 3-way

# LWafula 14082017
# hoteling's statistic method on glog transformed data: only on active treatments
# pick the 1st 10 features [for a dataset of 18 observations, 2 df lost becoz of means vecctors for each treatment in the pair, residual rank=16
# pick number of features sufficiently lower than 16 there4 to increase power of the MANOVA-based statistic
# this results are to be used to compare with the untranformed plot for the Hotelling's results and see if there's a better separation
# alpha estimates are available from the LWafula GLOG parameters estimation script in the MThesis/scripts folder
# the y-alpha dataset is also built in the LWafula GLOG parameters estimation script

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

datawd<-'E:/uhasselt/22yrsem/MThesis/data'
datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy'
datawdmatrixpvals<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/results-matricesH2_pvalues etc'

# 650 treatments, top 10 features
load('E:/uhasselt/22yrsem/MThesis/data/p1602xxPPS_means_preprocessed.Rdata')
datpOpt.feature.orig<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features_mRMR_optN")[1:10])]
rm(list=c('datp'))

# the transformed data [cell and well- level for lambda = 10.5]
lambv<-c(0.1,0.5,8.5,10.5)
for (i in 1:length(lambv)){
  #cell-level transformed data
  fileb1<-paste0("p1602xxPPS_means_preprocessedfeat_act_",lambv[i],'.RData')
  fileb1<-file.path(datawdtranslambda,fileb1)
  #lambda value 
  lambdaval=lambv[i]; lambdaval
  
  load(fileb1)
  datpopt.featurelambda<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(glogdatp, "features_mRMR_optN")[1:10])]
  rm(list=c('glogdatp','fileb1'))
  
  # well level
  fileb1<-paste0("aggregated_p1602xxPPS_means_preprocessedfeat_act_",lambv[i],'.RData')
  fileb1<-file.path(datawdtranslambda,fileb1)
  load(fileb1)
  welldatpopt.featurelambda<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(glogdatp, "features_mRMR_optN")[1:10])]
  rm(list=c('glogdatp','fileb1'))
  
  # for all, just keep those treatments that were active in both untransformed and transformed [cell and well level] data sets
  int.all<-Reduce(intersect, list(names(table(datpOpt.feature.orig$Treatment)),names(table(datpopt.featurelambda$Treatment)),
                                  names(table(welldatpopt.featurelambda$Treatment)))) # 628 for lambda= 0.5
  
  #datpOpt.feature.temp<-datpOpt.feature
  datpOpt.feature<- datpOpt.feature.orig[datpOpt.feature.orig$Treatment %in% int.all,]
  datpopt.featurelambda<- datpopt.featurelambda[datpopt.featurelambda$Treatment %in% int.all,]
  welldatpopt.featurelambda<-welldatpopt.featurelambda[welldatpopt.featurelambda$Treatment %in% int.all,]
  
  # ### ````````````````` untransformed for that specific lambda ``````````````````````
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
  
  results_H2_allTrts10fts_untransformed<-
    paste0('report/results-matricesH2_pvalues etc/cellwell/results_untransformed_lambda',lambdaval,'.csv')
  write.csv(results,results_H2_allTrts10fts_untransformed)
  
  # pvalues
  # into one line for plotting
  results_pvalue<-melt(results_pvalue)
  results_pvalue<-results_pvalue[complete.cases(results_pvalue[['value']]),]
  results_H2pvalues_allTrts10fts_untransformed<-
    paste0('report/results-matricesH2_pvalues etc/cellwell/pvalues_untransformed_lambda',lambdaval,'.csv')
  write.csv(results_pvalue,results_H2pvalues_allTrts10fts_untransformed)
  
  
  # ~~~~~~~~~~~~~~~~~~~~  cell level transformed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # distribution of the Hoteling's T2 statistic and its p-value
  treatkp<-names(table(datpopt.featurelambda$Treatment))
  results<-array(NA,dim=c(length(names(table(datpopt.featurelambda$Treatment))),length(names(table(datpopt.featurelambda$Treatment)))))
  colnames(results)<-paste(treatkp,sep = '')
  rownames(results)<-paste(treatkp,sep = '')
  
  #distribution of pvalues
  results_pvalue<-array(NA,dim=c(length(names(table(datpopt.featurelambda$Treatment))),length(names(table(datpopt.featurelambda$Treatment)))))
  colnames(results_pvalue)<-paste(treatkp,sep = '')
  rownames(results_pvalue)<-paste(treatkp,sep = '')
  
  # NOTE: thz takes time (2hrs), run without assigning the names if needed, otherwise use the saved results
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
  results<-melt(results)
  results<-results[complete.cases(results[['value']]),]
  
  results_H2_allTrts10fts_transformed<-
    paste0('report/results-matricesH2_pvalues etc/cellwell/results_celltransformed_lambda',lambdaval,'.csv')
  write.csv(results,results_H2_allTrts10fts_transformed)
  
  # pvalues
  results_pvalue<-melt(results_pvalue)
  results_pvalue<-results_pvalue[complete.cases(results_pvalue[['value']]),]
  results_H2pvalues_allTrts10fts_transformed<-
    paste0('report/results-matricesH2_pvalues etc/cellwell/pvalues_celltransformed_lambda',lambdaval,'.csv')
  write.csv(results_pvalue,results_H2pvalues_allTrts10fts_transformed)
  
  # ~~~~~~~~~~~~~~~~~~~~  well level transformed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # distribution of the Hoteling's T2 statistic and its p-value
  treatkp<-names(table(welldatpopt.featurelambda$Treatment))
  results<-array(NA,dim=c(length(names(table(welldatpopt.featurelambda$Treatment))),length(names(table(welldatpopt.featurelambda$Treatment)))))
  colnames(results)<-paste(treatkp,sep = '')
  rownames(results)<-paste(treatkp,sep = '')
  
  #distribution of pvalues
  results_pvalue<-array(NA,dim=c(length(names(table(welldatpopt.featurelambda$Treatment))),length(names(table(welldatpopt.featurelambda$Treatment)))))
  colnames(results_pvalue)<-paste(treatkp,sep = '')
  rownames(results_pvalue)<-paste(treatkp,sep = '')
  
  # NOTE: thz takes time (2hrs), run without assigning the names if needed, otherwise use the saved results
  for (j in 1:length(treatkp)) {
    for (i in 1:length(treatkp)) {
      if(j<i) { #you may as well do for (j!=i) though it is just a mirror-image
        fit_ij<-welldatpopt.featurelambda[welldatpopt.featurelambda$Treatment %in% treatkp[c(j,i)],]
        
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
  
  results_H2_allTrts10fts_transformed<-
    paste0('report/results-matricesH2_pvalues etc/cellwell/results_welltransformed_lambda',lambdaval,'.csv')
  write.csv(results,results_H2_allTrts10fts_transformed)
  
  # pvalues
  results_pvalue<-melt(results_pvalue)
  results_pvalue<-results_pvalue[complete.cases(results_pvalue[['value']]),]
  results_H2pvalues_allTrts10fts_transformed<-
    paste0('report/results-matricesH2_pvalues etc/cellwell/pvalues_welltransformed_lambda',lambdaval,'.csv')
  write.csv(results_pvalue,results_H2pvalues_allTrts10fts_transformed)
}


# 14-08-2017: ALL THE ABOVE ARE ALREADY SAVED- every lambda has its corresponding untranformed values and cell- and well-level results
# for treatments actiively called in all three sets
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

cellwell<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/results-matricesH2_pvalues etc/cellwell'
lambv<-c(0.1,0.5,8.5,10.5)
for (i in 1:length(lambv)){
  #                                         untransformed
  untransformed<-file.path(cellwell, paste0('results_untransformed_lambda',lambv[i],'.csv'))
  untransformed<-read.table(untransformed,sep=',',header = T)
  untransformed<-untransformed[,-1]
  myhist <- hist(untransformed$value)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(untransformed$value)
  mydensity$y <- mydensity$y * multiplier[1]
  
  # data.frame so we can add the rest
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'0'
  #                                            single-cell level
  cell<-file.path(cellwell, paste0('results_celltransformed_lambda',lambv[i],'.csv'))
  cell<-read.table(cell,sep=',',header = T)
  cell<-cell[,-1]
  
  #add the lines
  cellhisttemp <- hist(cell$value,plot=FALSE)
  multiplier <- cellhisttemp$counts / cellhisttemp$density
  celldensitytemp <- density(cell$value)
  celldensitytemp$y <- celldensitytemp$y * multiplier[1]
  celldensitytemp<-data.frame(cbind(y=celldensitytemp$y,x=celldensitytemp$x))
  celldensitytemp$gr<-'1'
  mydensity<-data.frame(rbind(mydensity,celldensitytemp))
  rm(celldensitytemp)
  
  #                                          well level 
  well<-file.path(cellwell, paste0('results_welltransformed_lambda',lambv[i],'.csv'))
  well<-read.table(well,sep=',',header = T)
  well<-well[,-1]
  
  #add the lines
  wellhisttemp <- hist(well$value,plot=FALSE)
  multiplier <- wellhisttemp$counts / wellhisttemp$density
  welldensitytemp <- density(well$value)
  welldensitytemp$y <- welldensitytemp$y * multiplier[1]
  welldensitytemp<-data.frame(cbind(y=welldensitytemp$y,x=welldensitytemp$x))
  welldensitytemp$gr<-'2'
  mydensity<-data.frame(rbind(mydensity,welldensitytemp))
  rm(welldensitytemp)
  
  #      plot distribution of respective Hotelling's T2
  mydensity$gr<-as.numeric(mydensity$gr)
  gr<-sort(unique(mydensity$gr))
  
  #plots
  hi<-file.path(cellwell,paste0('graphs/D2_matched_aggregated_hotellings_allactiveTRTS_lambda=',lambv[i],'.pdf'))
  pdf(hi)
  
  plot(mydensity[mydensity$gr==gr[1],]$y~mydensity[mydensity$gr==gr[1],]$x,type='l',lwd=1,lty=2,ylim=c(0,max(mydensity$y)),
       col='black',xlab=expression(paste('T'^2,' ','value')),ylab='Density',
       main=expression(paste('Distribution of',' ','T'^2)))
  lines(mydensity[mydensity$gr==1,]$y~mydensity[mydensity$gr==1,]$x,col='black',lwd=1,lty=1) #single cell level
  lines(mydensity[mydensity$gr==2,]$y~mydensity[mydensity$gr==2,]$x,col='red',lwd=1,lty=1) #welllevel
  legend('topright',legend = 
           c(expression(paste('Untransformed')),eval(bquote(expression(paste(lambda ~'='~.(lambv[i])~'[cell level]')))),
             eval(bquote(expression(paste(lambda ~'='~.(lambv[i])~'[well level]'))))),col=c('black','black','red'),bty = 'n',cex=0.85,lwd=c(1,1,1),lty=c(2,1,1))
  dev.off()
  
  #                                           logs
  #                                         untransformed
  myhistlog <- hist(log(untransformed$value))
  multiplier <- myhistlog$counts #/ myhistlog$density
  mydensitylog <- density(log(untransformed$value))
  mydensitylog$y <- mydensitylog$y * multiplier[1]
  
  # data.frame so we can add the rest
  mydensitylog<-data.frame(cbind(y=mydensitylog$y,x=mydensitylog$x))
  mydensitylog$gr<-'0'
  #                                            single-cell level
  
  cellhisttemplog <- hist(log(cell$value),plot=FALSE)
  multiplier <- cellhisttemplog$counts # / cellhisttemplog$density
  celldensitytemplog <- density(log(cell$value))
  celldensitytemplog$y <- celldensitytemplog$y * multiplier[1]
  celldensitytemplog<-data.frame(cbind(y=celldensitytemplog$y,x=celldensitytemplog$x))
  celldensitytemplog$gr<-'1'
  mydensitylog<-data.frame(rbind(mydensitylog,celldensitytemplog))
  rm(celldensitytemplog)
  
  #                                          well level 
  wellhisttemplog <- hist(log(well$value),plot=FALSE)
  multiplier <- wellhisttemplog$counts # / wellhisttemplog$density
  welldensitytemplog <- density(log(well$value))
  welldensitytemplog$y <- welldensitytemplog$y * multiplier[1]
  welldensitytemplog<-data.frame(cbind(y=welldensitytemplog$y,x=welldensitytemplog$x))
  welldensitytemplog$gr<-'2'
  mydensitylog<-data.frame(rbind(mydensitylog,welldensitytemplog))
  rm(welldensitytemplog)
  
  #                                                 plots
  mydensitylog$gr<-as.numeric(mydensitylog$gr)
  gr<-sort(unique(mydensity$gr))
  
  # plot
  hi<-file.path(cellwell,paste0('graphs/D2_LOG_matched_aggregated_hotellings_allactiveTRTS_lambda=',lambv[i],'.pdf'))
  pdf(hi)
  
  plot(mydensitylog[mydensitylog$gr==gr[1],]$y~mydensitylog[mydensitylog$gr==gr[1],]$x,type='l',lwd=1,lty=1,ylim=c(0,max(mydensitylog$y)),
       xlim=c(0,14),col='black',xlab=expression(paste('log(',' ','T'^2,' ','value',')')),ylab='Density',
       main=expression(paste('Distribution of',' ','log(T'^2,')')))
  lines(mydensitylog[mydensitylog$gr==1,]$y~mydensitylog[mydensitylog$gr==1,]$x,col='blue',lwd=1,lty=1) #single cell level
  lines(mydensitylog[mydensitylog$gr==2,]$y~mydensitylog[mydensitylog$gr==2,]$x,col='red',lwd=1,lty=1) #well level
  abline(v=median(mydensitylog[mydensitylog$gr==0,]$x),lwd=1,lty=2, col='black') #median for untransformed
  abline(v=median(mydensitylog[mydensitylog$gr==1,]$x),lwd=1,lty=2, col='blue') #median for single cell level
  abline(v=median(mydensitylog[mydensitylog$gr==2,]$x),lwd=1,lty=2, col='red') #median for well level
  
  legend('topleft',legend = 
           c(expression(paste('Untransformed')),eval(bquote(expression(paste(lambda ~'='~.(lambv[i])~'[cell level]')))),
             eval(bquote(expression(paste(lambda ~'='~.(lambv[i])~'[well level]')))),
             expression(paste('median untransformed')),expression(paste('median log(','T'^2,')',' ','transformed-single cell')),
             expression(paste('median log(','T'^2,')',' ','transformed-well level'))),col=c('black','blue','red','black','blue','red'),
         bty = 'n',cex=0.85,lwd=c(1,1,1,1,1,1), lty=c(1,1,1,2,2,2))
  
  dev.off()
  #                     Boxplots
  mydensitylog$gr<-factor(mydensitylog$gr,levels = c(0,1,2),labels = c('untransformed','cell-level','well-level'))
  
  hi<-file.path(cellwell,paste0('graphs/D2_LOG_matched_aggregated_Boxplothotellings_allactiveTRTS_lambda=',lambv[i],'.pdf'))
  pdf(hi)
  par(xpd = T, mar = par()$mar + c(0,1,0,1))
  plot(mydensitylog$x~as.factor(mydensitylog$gr),ylab=expression(paste('log(',' ','T'^2,' ','value',')')),xlab='')
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
}









