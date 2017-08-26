
# DATASET 2
### MATCHED: comparison for Hotelling's T2 for only treatments that were active pre- and post-transformation [using lambda= 10.5]


# LWafula 13082017
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

datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy/for H2'
datawdmatrixpvals<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/results-matricesH2_pvalues etc'

# 650 treatments, top 10 features
load('E:/uhasselt/22yrsem/MThesis/data/p1602xxPPS_means_preprocessed.Rdata')
datpOpt.feature.orig<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features_mRMR_optN")[1:10])]
rm(list=c('datp'))

# the transformed data
file.names <- dir(datawdtranslambda, pattern ="p1602xxPPS_means_preprocessedfeat_act_")
file.names<-file.names[!grepl('(aggregated_)',file.names, ig = TRUE)]

# for all the lambdas, update with matched Hotelling's T2
for (k in 1:length(file.names)){
  fileb1<-file.path(datawdtranslambda,file.names[k])
  #lambda value 
  lambdaval=gsub("p1602xxPPS_means_preprocessedfeat_act_","",file.names[k]); lambdaval<-gsub(".RData","",lambdaval); lambdaval
  
  load(fileb1)
  datpopt.featurelambda<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(glogdatp, "features_mRMR_optN")[1:10])]
  rm(list=c('glogdatp','fileb1'))
  
  # for both, just keep those treatments that were active in both untransformed and transformed data sets
  intersect(names(table(datpOpt.feature.orig$Treatment)),names(table(datpopt.featurelambda$Treatment)) ) # lambda specific, e.g lambda=10.5,they r 683
  
  #datpOpt.feature.temp<-datpOpt.feature
  datpOpt.feature<- datpOpt.feature.orig[datpOpt.feature.orig$Treatment %in% 
                                           intersect(names(table(datpOpt.feature.orig$Treatment)),names(table(datpopt.featurelambda$Treatment))),]
  datpopt.featurelambda<- datpopt.featurelambda[datpopt.featurelambda$Treatment %in% 
                                                  intersect(names(table(datpOpt.feature.orig$Treatment)),names(table(datpopt.featurelambda$Treatment))),]
  
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
    paste0('report/results-matricesH2_pvalues etc/matchedactivetreatments/results_untransformed_lambda',lambdaval,'.csv')
  write.csv(results,results_H2_allTrts10fts_untransformed)
  
  # pvalues
  # into one line for plotting
  results_pvalue<-melt(results_pvalue)
  results_pvalue<-results_pvalue[complete.cases(results_pvalue[['value']]),]
  results_H2pvalues_allTrts10fts_untransformed<-
    paste0('report/results-matricesH2_pvalues etc/matchedactivetreatments/pvalues_untransformed_lambda',lambdaval,'.csv')
  write.csv(results_pvalue,results_H2pvalues_allTrts10fts_untransformed)
  
  
  # ~~~~~~~~~~~~~~~~~~~~ transformed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    paste0('report/results-matricesH2_pvalues etc/matchedactivetreatments/results_transformed_lambda',lambdaval,'.csv')
  write.csv(results,results_H2_allTrts10fts_transformed)
  
  # pvalues
  results_pvalue<-melt(results_pvalue)
  results_pvalue<-results_pvalue[complete.cases(results_pvalue[['value']]),]
  results_H2pvalues_allTrts10fts_transformed<-
    paste0('report/results-matricesH2_pvalues etc/matchedactivetreatments/pvalues_transformed_lambda',lambdaval,'.csv')
  write.csv(results_pvalue,results_H2pvalues_allTrts10fts_transformed)
  
  #
  print(lambdaval)
}



# distribution of T2 for all lambdas + lOGS
# each lambda has its own version of untransformed results for the T2 since matched active treatments vary across transformations and there4
# commonly called treatments btn untransformed and lambda-specific transformations differ

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
matched<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/results-matricesH2_pvalues etc/matchedactivetreatments'
lambv<-c(0.1,seq(0.5,25,by=0.5))
data<-c()
for(i in 1:length(lambv)){
  #                                          lambda[i]-specific T2 for untransformed
  untransformed<-file.path(matched,paste0('results_untransformed_lambda',lambv[i],'.csv'))
  untransformed<-read.table(untransformed,sep=',',header = TRUE)
  untransformed<-untransformed[,-1]
  untransformed$lambda<-lambv[i]
  untransformed$gr<-'0'
  untransformed<-untransformed[,c('value','lambda','gr')]
  
  #                                    transformed
  transformed<-file.path(matched,paste0('results_transformed_lambda',lambv[i],'.csv'))
  transformed<-read.table(transformed,sep=',',header = TRUE)
  transformed<-transformed[,-1]
  transformed$lambda<-lambv[i]
  transformed$gr<-'1'
  transformed<-transformed[,c('value','lambda','gr')]
  
  data<-rbind(data,untransformed,transformed)
}
data$gr<-as.numeric(data$gr)
datasave<-'report/results-matricesH2_pvalues etc/matchedactivetreatments/data2Boxplots.RData'
save(data,file=datasave)

#http://www.r-graph-gallery.com/9-ordered-boxplot/
# reorder the lambdas from the most sensible to the most resistant. (mixing low and high treatments for the calculations)
new_order <- with(data, reorder(lambda , value, mean , na.rm=T))

# Then I make the boxplot, asking to use the 2 factors : variety (in the good order) AND treatment :
pdf<-'report/dataView/manova plots/individual lambdas/D2_matched_loghotellings_allactiveTRTS_Boxplot_all_lambdas.pdf'
pdf('pdf')
par(mar=c(4.5,4.5,3,1))
myplot=boxplot(log(value) ~ gr*lambda , data=data  , boxwex=0.75 , ylab=expression(paste('log(',' ','T'^2,' ','value',')')),
               main="",xlab=expression(paste(lambda,' ','value')) , col=c("black" , "red") ,  xaxt="n",pch=19,cex=0.25,
               ylim=c(-1,max(log(data$value))+1))

# To add the label of x axis
strsplit(myplot$names , '\\.')->lambdas
my_lambdas<-c()
for(i in 1:length(lambdas)){
  if(length(lambdas[[i]])==2){
    my_lambdas[i]<-lambdas[[i]][2]
  }
  if(length(lambdas[[i]])==3){
    my_lambdas[i]<-paste0(lambdas[[i]][2],'.',lambdas[[i]][3])
  }
}
my_lambdas=my_lambdas[seq(1 , length(my_lambdas) , 2)]; my_lambdas<-sort(as.numeric(my_lambdas))
axis(1, at = seq(1.5 , 2*length(lambv) , 2), labels = my_lambdas , tick=TRUE , cex=0.3)
for(i in seq(0.35,2*length(lambv),4)){ abline(v=i,lty=1, col="black")}
# Add a legend
legend("bottomright", legend = c("untransformed", "transformed"), col=c("black" , "red"),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  ncol = 2, inset = c(0.1, 0.1))
dev.off()
rm(untransformed,transformed)

# distribution of individual lambdas vs respective untransformed T2
for(i in 1:length(lambv)){
  
  #                                         untransformed
  untransformed<-data[data$gr==0 & data$lambda==lambv[i] ,]
  myhist <- hist(untransformed$value)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(untransformed$value)
  mydensity$y <- mydensity$y * multiplier[1]
  
  # data.frame so we can add the rest
  mydensity<-data.frame(cbind(y=mydensity$y,x=mydensity$x))
  mydensity$gr<-'0'
  #                                            single-cell level
  transformed<-data[data$gr==1 & data$lambda==lambv[i] ,]
  
  #add the lines
  transformedhisttemp <- hist(transformed$value,plot=FALSE)
  multiplier <- transformedhisttemp$counts / transformedhisttemp$density
  transformeddensitytemp <- density(transformed$value)
  transformeddensitytemp$y <- transformeddensitytemp$y * multiplier[1]
  transformeddensitytemp<-data.frame(cbind(y=transformeddensitytemp$y,x=transformeddensitytemp$x))
  transformeddensitytemp$gr<-'1'
  mydensity<-data.frame(rbind(mydensity,transformeddensitytemp))
  rm(transformeddensitytemp)
  
  #      plot distribution of respective Hotelling's T2
  mydensity$gr<-as.numeric(mydensity$gr)
  gr<-sort(unique(mydensity$gr))
  
  #plots
  hi<-file.path(matched,paste0('graphs/D2_matched_hotellings_allactiveTRTS_lambda=',lambv[i],'.pdf'))
  pdf(hi)
  
  plot(mydensity[mydensity$gr==gr[1],]$y~mydensity[mydensity$gr==gr[1],]$x,type='l',lwd=2,lty=1,ylim=c(0,max(mydensity$y)),
       col='black',xlab=expression(paste('T'^2,' ','value')),ylab='Density',
       main=expression(paste('Distribution of',' ','T'^2)))
  lines(mydensity[mydensity$gr==1,]$y~mydensity[mydensity$gr==1,]$x,col='blue',lwd=2,lty=1) #single cell level
  legend('topright',legend = 
           c(expression(paste('Untransformed')),eval(bquote(expression(paste(lambda ~'='~.(lambv[i])))))),
         col=c('black','blue'),bty = 'n',cex=0.85,lwd=c(2,2),lty=c(1,1))
  dev.off()
  
  #clipped for better visuals
  hi<-file.path(matched,paste0('graphs/D2_clipped_matched_hotellings_allactiveTRTS_lambda=',lambv[i],'.pdf'))
  pdf(hi)
  
  plot(mydensity[mydensity$gr==gr[1],]$y~mydensity[mydensity$gr==gr[1],]$x,type='l',lwd=2,lty=1,ylim=c(0,max(mydensity$y)),
       col='black',xlab=expression(paste('T'^2,' ','value')),ylab='Density',
       main=expression(paste('Distribution of',' ','T'^2)),xlim=c(min(mydensity$x),100000))
  lines(mydensity[mydensity$gr==1,]$y~mydensity[mydensity$gr==1,]$x,col='blue',lwd=2,lty=1) #single cell level
  legend('topright',legend = 
           c(expression(paste('Untransformed')),eval(bquote(expression(paste(lambda ~'='~.(lambv[i])))))),
         col=c('black','blue'),bty = 'n',cex=0.85,lwd=c(2,2),lty=c(1,1))
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
  #                                            single-cell transformed
  
  transformedhisttemplog <- hist(log(transformed$value),plot=FALSE)
  multiplier <- transformedhisttemplog$counts # / transformedhisttemplog$density
  transformeddensitytemplog <- density(log(transformed$value))
  transformeddensitytemplog$y <- transformeddensitytemplog$y * multiplier[1]
  transformeddensitytemplog<-data.frame(cbind(y=transformeddensitytemplog$y,x=transformeddensitytemplog$x))
  transformeddensitytemplog$gr<-'1'
  mydensitylog<-data.frame(rbind(mydensitylog,transformeddensitytemplog))
  rm(transformeddensitytemplog)
  
  #                                                 plots
  mydensitylog$gr<-as.numeric(mydensitylog$gr)
  gr<-sort(unique(mydensity$gr))
  
  # plot
  hi<-file.path(matched,paste0('graphs/D2_LOG_matched__hotellings_allactiveTRTS_lambda=',lambv[i],'.pdf'))
  pdf(hi)
  
  plot(mydensitylog[mydensitylog$gr==gr[1],]$y~mydensitylog[mydensitylog$gr==gr[1],]$x,type='l',lwd=2,lty=1,ylim=c(0,max(mydensitylog$y)+0.05),xlim=c(0,14),
       col='black',xlab=expression(paste('log(',' ','T'^2,' ','value',')')),ylab='Density',
       main=expression(paste('Distribution of',' ','log(T'^2,')')))
  lines(mydensitylog[mydensitylog$gr==1,]$y~mydensitylog[mydensitylog$gr==1,]$x,col='blue',lwd=2,lty=1) #single cell level
  abline(v=median(mydensitylog[mydensitylog$gr==0,]$x),lwd=2,lty=2, col='black') #median for untransformed
  abline(v=median(mydensitylog[mydensitylog$gr==1,]$x),lwd=2,lty=2, col='blue') #median for single cell level
  
  legend('topleft',legend = 
           c(expression(paste('Untransformed')),eval(bquote(expression(paste(lambda ~'='~.(lambv[i])~'[cell level]')))),
             expression(paste('median untransformed')),expression(paste('median transformed'))),
         col=c('black','blue','black','blue'),
         bty = 'n',cex=0.85,lwd=c(2,2,2,2), lty=c(1,1,2,2))
  
  dev.off()
  #                     Boxplots
  mydensitylog$gr<-factor(mydensitylog$gr,levels = c(0,1),labels = c('untransformed','transformed'))
  
  hi<-file.path(matched,paste0('graphs/D2_LOG_matched_Boxplothotellings_allactiveTRTS_lambda=',lambv[i],'.pdf'))
  pdf(hi)
  par(xpd = T, mar = par()$mar + c(0,1,0,1))
  plot(mydensitylog$x~as.factor(mydensitylog$gr),ylab=expression(paste('log(',' ','T'^2,' ','value',')')),xlab='')
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
}

rm(untransformed,transformed)


















