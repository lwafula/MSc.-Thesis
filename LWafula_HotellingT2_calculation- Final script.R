#DATASET2

# LWafula 27072017
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
datawdtranslambda<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda - Copy/for H2'
datawdmatrixpvals<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/results-matricesH2_pvalues etc'

# for the untransformed data
fn_out0          <- 'p1602xxPPS_means_preprocessedfeat_act.Rdata'
fn_out0          <- file.path(datawd, fn_out0)
load(fn_out0)

# 650 treatments, top 10 features
load('E:/uhasselt/22yrsem/MThesis/data/p1602xxPPS_means_preprocessed.Rdata')
datpOpt.feature<-datp[datp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(datp, "features_mRMR_optN")[1:10])]
rm(list=c('datp','fn_out0','datpfeatact'))

# distribution of the Hoteling's T2 statistic and its p-value
treatkp<-names(table(datpOpt.feature$Treatment))
results<-array(NA,dim=c(length(names(table(datpOpt.feature$Treatment))),length(names(table(datpOpt.feature$Treatment)))))

colnames(results)<-paste(names(table(datpOpt.feature$Treatment)),sep = '')
rownames(results)<-paste(names(table(datpOpt.feature$Treatment)),sep = '')

#distribution of pvalues
results_pvalue<-array(NA,dim=c(length(names(table(datpOpt.feature$Treatment))),length(names(table(datpOpt.feature$Treatment)))))

colnames(results_pvalue)<-paste(names(table(datpOpt.feature$Treatment)),sep = '')
rownames(results_pvalue)<-paste(names(table(datpOpt.feature$Treatment)),sep = '')

# NOTE: thz takes time (2hrs), run without assigning the names if needed, otherwise use the saved results
for (j in 1:length(treatkp)) {
  for (i in 1:length(treatkp)) {
    if(j<i) { #you may as well do for (j!=i) though it is just a mirror-image
      fit_ij<-datpOpt.feature[datpOpt.feature$Treatment %in% treatkp[c(j,i)],]
      treatkp<-names(table(datpOpt.feature$Treatment))
      
      #calculating the Hotelling and pvalue
      (fit<-hotelling.test(.~Treatment, data = fit_ij[,c(1,2:11)]))
      #T2
      results[j,i]<-fit$stats$statistic
      #pval
      results_pvalue[j,i]<-fit$pval
      print(fit$stats$statistic)
    }
    else {print(paste('i=',i,'j=',j,'already done/no need to do'))}
  }
}

# How to plot a table of numbers such that the values are represented by color?
pdf('report/dataView/manova plots/D2_imageT2_alltreats_untransformed.pdf')
image(results,main='Untransformed data')
dev.off()

#save the matrices b4 proceeding
results_H2_allTrts10fts_untransformed<-'report/results-matricesH2_pvalues etc/July2017results_H2_allTrts10fts_untransformed.csv'
results<-melt(results); results<-results[complete.cases(results$value),]
write.csv(results,results_H2_allTrts10fts_untransformed)
# #
results_H2pvalues_allTrts10fts_untransformed<-'report/results-matricesH2_pvalues etc/July2017results_H2pvalues_allTrts10fts_untransformed.csv'
results_pvalue<-melt(results_pvalue); results_pvalue<-results_pvalue[complete.cases(results_pvalue$value),]
write.csv(results_pvalue,results_H2pvalues_allTrts10fts_untransformed)

# hoteling's dist. plots

# which treatment combination has the largest difference? Check the best performing lambda, pick out trt combinations wth large difference & plot pre-post profiles
# View(results[results$value==summary(results$value)[6],]) 
# [cmpd892 @ 9 uM in Hepa-ER-Mito & cmpd94 @ 9 uM in Hepa-ER-Mito with T2=242201.3]
# no diff really: cmpd417 @ 1 uM in Hepa-ER-Mito & cmpd417 @ 3 uM in Hepa-ER-Mito with T2=2.42051

myhist <- hist(results$value)
multiplier <- myhist$counts / myhist$density
mydensity <- density(results$value)
mydensity$y <- mydensity$y * multiplier[1]

pdf('report/dataView/manova plots/D2_hotellings_allactiveTRTS_untransformed.pdf')
plot(myhist,xlab=expression(paste('T'^2,' ','value')),ylab='counts',xlim=c(0,350000),main=paste('calculated Hotelling',"'",'s statistic- Untransformed data'),
     ylim=c(0,1000000))
lines(mydensity,lwd=2)
dev.off()

#add here plots for other transformations for the report and H2 comparisons btwn untransformed and transformed
plot(mydensity,lwd=2,col='black',xlab=expression(paste('T'^2,' ','value')),main='',xlim=c(0,350000))

# pvalues

myhist_pvalue <- hist(results_pvalue$value)
multiplier_pvalue <- myhist_pvalue$counts / myhist_pvalue$density
mydensity_pvalue <- density(results_pvalue$value)
mydensity_pvalue$y <- mydensity_pvalue$y * multiplier_pvalue[1]

plot(myhist_pvalue,xlab=expression(paste('P-',' ','value')),ylab='counts',xlim=c(0,0.25))

# add results for lambdas [experimented using lambda=2.5] then loop through all available lambdas and make a plot
# note that some datasets were initially saved wrongly as p1602.. script corrected as of now [26-07-2017]
file.names <- dir(datawdtranslambda, pattern ="p1602xxPPS_means_preprocessedfeat_act_")
for (i in 1:length(file.names)){
  fileb1<-file.path(datawdtranslambda,file.names[i])
  load(fileb1)
  datpopt.featurelambda<-glogdatp[glogdatp[['Active_FractionOfRepl']]>=0.5 ,c("Treatment", attr(glogdatp, "features_mRMR_optN")[1:10])]
  rm(list=c('glogdatp','fileb1'))

  # distribution of the Hoteling's T2 statistic and its p-value
  treatkp<-names(table(datpopt.featurelambda$Treatment))
  #lambda value
  lambdaval=   gsub("p1602xxPPS_means_preprocessedfeat_act_","",file.names[i]); lambdaval<-gsub(".RData","",lambdaval); lambdaval
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
  #image plot
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
  
  #save the matrices b4 proceeding
  resultsname<-paste0('results',lambdaval,'.csv')
  resultsname<-file.path(datawdmatrixpvals,resultsname)
  write.csv(results,resultsname)
  
  results_pvalue<-melt(results_pvalue)
  results_pvalue<-results_pvalue[complete.cases(results_pvalue[['value']]),]
  
  resultsnamepvalues<-paste0('results_pvalue',lambdaval,'.csv')
  resultsnamepvalues<-file.path(datawdmatrixpvals,resultsnamepvalues)
  write.csv(results_pvalue,resultsnamepvalues)
  
  myhist <- hist(results$value)
  multiplier <- myhist$counts / myhist$density
  mydensity <- density(results$value)
  mydensity$y <- mydensity$y * multiplier[1]
  
  h2<-paste0('report/dataView/manova plots/D2_hotellings_allactiveTRTS_lambda=',lambdaval,'.pdf')
  
  pdf(h2)
  plot(myhist,xlab=expression(paste('T'^2,' ','value')),ylab='counts',xlim=c(0,summary(results$value)[6]),
       main=paste('calculated Hotelling',"'",'s statistic- transformed data'),ylim=c(0,1000000))
  lines(mydensity,lwd=2)
  dev.off()
  
}



# 28-07-2017: ALL THE ABOVE ARE ALREADY SAVED
# Hotelling's distribution

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
  
  # into one line for plotting; already molten
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

#plots
hi<-paste0('report/dataView/manova plots/individual lambdas/D2_hotellings_allactiveTRTS_lambda=','all lambdas','.pdf')
pdf(hi)

plot(mydensity[mydensity$gr==gr[1],]$y~mydensity[mydensity$gr==gr[1],]$x,type='l',lwd=1,lty=2,ylim=c(0,max(mydensity$y)),
     col='black',xlab=expression(paste('T'^2,' ','value')),ylab='Density',
     main=expression(paste('Distribution of',' ','T'^2,' ','statistic with changing',' ',lambda)))

for(j in 1:length(gr)){
  if(gr[j]<1){
    lines(mydensity[mydensity$gr==gr[j],]$y~mydensity[mydensity$gr==gr[j],]$x,col='black',lwd=1,lty=2)
    legend('top',legend = expression(paste(lambda,' ','= 0 [Untransformed]')),col='black',lty=2,bty = 'n',cex=0.5,lwd=1) 
  }
  if(gr[j]>0 & gr[j]<=10){
    lines(mydensity[mydensity$gr==gr[j],]$y~mydensity[mydensity$gr==gr[j],]$x,col=j,lwd=2)
    legend(100000,((52e+05)*0.071*gr[j]+5e+05),legend = gr[j],col=unique(j),lty=1,bty = 'n',cex=0.5,lwd=2)
  }
   if(gr[j]>10 & gr[j]<=20){
     lines(mydensity[mydensity$gr==gr[j],]$y~mydensity[mydensity$gr==gr[j],]$x,col=j,lwd=2)
     legend(150000,((52e+05)*0.071*gr[j]-32e+05),legend = gr[j],col=unique(j),lty=1,bty = 'n',cex=0.5,lwd=2)
   }
   if(gr[j]>20){
     lines(mydensity[mydensity$gr==gr[j],]$y~mydensity[mydensity$gr==gr[j],]$x,col=j,lwd=2)
     legend(200000,((52e+05)*0.071*gr[j]-70e+05),legend = gr[j],col=unique(j),lty=1,bty = 'n',cex=0.5,lwd=2)
   }
   
}
dev.off()

############ logs for the above ############

# logs might be better

myhistlog <- hist(log(h2alltrst10features$value), density = NULL)
multiplierlog <- myhistlog$counts
mydensitylog <- density(log(h2alltrst10features$value))
mydensitylog$y <- mydensitylog$y * multiplierlog[1]

mydensitylog<-data.frame(cbind(y=mydensitylog$y,x=mydensitylog$x))
mydensitylog$gr<-'0'

# add transform distribution lines for Hotelling's
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
  h2temp<-h2temp[complete.cases(h2temp[['value']]),]
  
  #add the lines
  myhisttemp <- hist(log(h2temp$value),plot=FALSE)
  multiplier <- myhisttemp$counts
  mydensitytemp <- density(log(h2temp$value))
  mydensitytemp$y <- mydensitytemp$y * multiplier[1]
  
  mydensitytemp<-data.frame(cbind(y=mydensitytemp$y,x=mydensitytemp$x))
  mydensitytemp$gr<-lambv[i]
  mydensitylog<-data.frame(rbind(mydensitylog,mydensitytemp))
}
mydensitylog$gr<-as.numeric(mydensitylog$gr)
gr<-sort(unique(mydensitylog$gr))

# plot
  # #add here plots for other transformations for the report and H2 comparisons btwn untransformed and transformed
  hi<-paste0('report/dataView/manova plots/individual lambdas/logs/D2_LOG_hotellings_allactiveTRTS_lambda=','all lambdas','.pdf')
  pdf(hi)
  
  plot(mydensitylog[mydensitylog$gr==gr[1],]$y~mydensitylog[mydensitylog$gr==gr[1],]$x,type='l',lwd=1,lty=2,ylim=c(0,2.5),xlim=c(0,14),
       col='black',xlab=expression(paste('log(',' ','T'^2,' ','value',')')),ylab='Density',
       main=expression(paste('Distribution of',' ','log(T'^2,')',' ','statistic with changing',' ',lambda)))
  
  for(j in 1:length(gr)){
    if(gr[j]<1){
      lines(mydensitylog[mydensitylog$gr==gr[j],]$y~mydensitylog[mydensitylog$gr==gr[j],]$x,col='black',lwd=1,lty=2)
      legend('top',legend = expression(paste(lambda,' ','= 0 [Untransformed]')),col='black',bty = 'n',cex=0.5,lwd=1,lty=2) 
    }
    if(gr[j]>0 & gr[j]<=10){
      lines(mydensitylog[mydensitylog$gr==gr[j],]$y~mydensitylog[mydensitylog$gr==gr[j],]$x,col=j,lwd=2)
      legend(0,(gr[j]*0.2+0.2),legend = gr[j],col=unique(j),lty=1,bty = 'n',cex=0.5,lwd=2)
    }
    if(gr[j]>10 & gr[j]<=22.5){
      lines(mydensitylog[mydensitylog$gr==gr[j],]$y~mydensitylog[mydensitylog$gr==gr[j],]$x,col=j,lwd=2)
      legend(1.5,(gr[j]*0.15-1.25),legend = gr[j],col=unique(j),lty=1,bty = 'n',cex=0.5,lwd=2)
    }
    if (gr[j]>22.5){
      lines(mydensitylog[mydensitylog$gr==gr[j],]$y~mydensitylog[mydensitylog$gr==gr[j],]$x,col=j,lwd=2)
      legend(4,(gr[j]*0.2-4),legend = gr[j],col=unique(j),lty=1,bty = 'n',cex=0.5,lwd=2)
    }
  }
  dev.off()

########### end of logs for all ################

# Individuals
file.names <- dir(datawdmatrixpvals, pattern ="results")
#keep only Hotelling's results for now
file.names<-file.names[!grepl('^(results_)',file.names, ig = TRUE)]
file.names<-file.names[!grepl('^(results_)',file.names, ig = TRUE)]; file.names<-file.names[!grepl('^(July2017)',file.names, ig = TRUE)]

lambv<-c()

for (i in 1:length(file.names)){
  h2temp<-file.path(datawdmatrixpvals,file.names[i])
  lambv[i]<-gsub('results','',file.names[i]); lambv[i]<-gsub('.csv','',lambv[i]); lambv[i]
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  
  # into one line for plotting;   h2temp<-melt(h2temp)
  h2temp<-h2temp[complete.cases(h2temp[['value']]),]
  
  #add the lines
  myhisttemp <- hist(h2temp$value,plot=FALSE)
  multiplier <- myhisttemp$counts / myhisttemp$density
  mydensitytemp <- density(h2temp$value)
  mydensitytemp$y <- mydensitytemp$y * multiplier[1]
  
  mydensitytemp<-data.frame(cbind(y=mydensitytemp$y,x=mydensitytemp$x))
  mydensitytemp$gr<-lambv[i]
  
  #plots for each lambda
  hi<-paste0('report/dataView/manova plots/individual lambdas/D2_hotellings_allactiveTRTS_lambda=',lambv[i],'.pdf')
  pdf(hi)
  #look for a way of avoiding thz plot inside here!
  plot(mydensity[mydensity$gr==0,]$y~mydensity[mydensity$gr==0,]$x,type='l',lwd=1,lty=2,
       ylim=c(0,max(mydensity[mydensity$gr==0,]$y,mydensitytemp$y)),
       col='black',xlab=expression(paste('T'^2,' ','value')),ylab='Density',
       main=expression(paste('Distribution of',' ','T'^2,' ','statistic with changing',' ',lambda)))
  
  lines(mydensitytemp[mydensitytemp$gr==lambv[i],]$y~mydensitytemp[mydensitytemp$gr==lambv[i],]$x,col='blue',lwd=2)
  legend('topright',legend = c('untransformed',eval(bquote(expression(lambda~'='  ~.(lambv[i]))))),
         col=c('black','blue'),lty=c(2,1),bty = 'n',cex=1.5)
  dev.off()
}

# Individuals logs
file.names <- dir(datawdmatrixpvals, pattern ="results")
#keep only Hotelling's results for now
file.names<-file.names[!grepl('^(results_)',file.names, ig = TRUE)]
file.names<-file.names[!grepl('^(results_)',file.names, ig = TRUE)]; file.names<-file.names[!grepl('^(July2017)',file.names, ig = TRUE)]

lambv<-c()

for (i in 1:length(file.names)){
  h2temp<-file.path(datawdmatrixpvals,file.names[i])
  lambv[i]<-gsub('results','',file.names[i]); lambv[i]<-gsub('.csv','',lambv[i]); lambv[i]
  h2temp<-read.table(h2temp,header = T,sep=',')
  h2temp<-h2temp[,-1]
  
  # into one line for plotting -   h2temp<-melt(h2temp) [already done]
  h2temp<-h2temp[complete.cases(h2temp[['value']]),]
  
  #add the lines
  myhisttemp <- hist(log(h2temp$value),plot=FALSE)
  multiplier <- myhisttemp$counts
  mydensitytemp <- density(log(h2temp$value))
  mydensitytemp$y <- mydensitytemp$y * multiplier[1]
  
  mydensitytemp<-data.frame(cbind(y=mydensitytemp$y,x=mydensitytemp$x))
  mydensitytemp$gr<-lambv[i]
  
  #plots for each lambda
  hi<-paste0('report/dataView/manova plots/individual lambdas/logs/D2_LOGS_hotellings_allactiveTRTS_lambda=',lambv[i],'.pdf')
  pdf(hi)
  #look for a way of avoiding thz plot inside here!
  plot(mydensitylog[mydensitylog$gr==0,]$y~mydensitylog[mydensitylog$gr==0,]$x,type='l',lwd=1,lty=1,
       ylim=c(0,max(mydensitylog[mydensitylog$gr==0,]$y,mydensitytemp$y)),
       col='black',xlab=expression(paste('T'^2,' ','value')),ylab='Density',
       main=expression(paste('Distribution of',' ','T'^2,' ','statistic with changing',' ',lambda)))
  abline(v=mean(mydensitylog[mydensitylog$gr==0,]$x),lwd=1,lty=2, col='black')
  lines(mydensitytemp[mydensitytemp$gr==lambv[i],]$y~mydensitytemp[mydensitytemp$gr==lambv[i],]$x,col='red',lwd=2)
  abline(v=mean(mydensitytemp[mydensitytemp$gr==lambv[i],]$x),lwd=1,lty=2, col='red')
  legend('topleft',
         legend = c('untransformed',eval(bquote(expression(lambda~'='  ~.(lambv[i])))), expression(paste('mean log(','T'^2,')',' ','untransformed')),expression(paste('mean log(','T'^2,')',' ','transformed'))),
         col=c('black','red','black','red'),lty=c(1,1,2,2),bty = 'n',cex=1.5)
  dev.off()
}
## end of individual logs


# P-VALUES

results_H2pvalues_allTrts10fts_untransformed<-'July2017results_H2pvalues_allTrts10fts_untransformed.csv'
results_H2pvalues_allTrts10fts_untransformed<-file.path(datawdmatrixpvals,results_H2pvalues_allTrts10fts_untransformed)
results_H2pvalues_allTrts10fts_untransformed<-read.table(results_H2pvalues_allTrts10fts_untransformed,header=T,sep=',')
results_H2pvalues_allTrts10fts_untransformed<-results_H2pvalues_allTrts10fts_untransformed[,-1]

#results_H2pvalues_allTrts10fts_untransformed<-melt(results_H2pvalues_allTrts10fts_untransformed)
results_H2pvalues_allTrts10fts_untransformed<-results_H2pvalues_allTrts10fts_untransformed[complete.cases(results_H2pvalues_allTrts10fts_untransformed$value),]
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
  
  # into one line for plotting - [h2temp<-melt(h2temp)]
  h2temp<-h2temp[complete.cases(h2temp[['value']]),]
  h2temp$grp<-paste('lambda',lambv[i])
  h2temp<-h2temp[,-c(1,2)]
  
  #merge for plotting and differentiating
  merge<-rbind(results_H2pvalues_allTrts10fts_untransformed,h2temp)
  #scatter plots for each lambda against the untransformed
  hi<-paste0('report/dataView/manova plots/individual lambdas/D2_pvalue_hotellings_allactiveTRTS_lambda=',lambv[i],'.pdf')
  pdf(hi)
  par(xpd = T, mar = par()$mar + c(0,0,0,7))
  plot(merge[['value']],type='n',xlab='index',ylab='p-value',
  xlim=c(0,max(length(merge[merge[['grp']]=='untransformed',][['value']]),length(merge[merge[['grp']]==paste('lambda',lambv[i]),][['value']]))))
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
  legend(max(length(merge[merge[['grp']]=='untransformed',][['value']]),length(merge[merge[['grp']]==paste('lambda',lambv[i]),][['value']]))+
           0.05*(max(length(merge[merge[['grp']]=='untransformed',][['value']]),length(merge[merge[['grp']]==paste('lambda',lambv[i]),][['value']]))),0.8, legend = c(
    paste('untransformed =',unt_pval5),paste('transformed =',t_pval5),'0.05 significance level'),pch=c(19,19,NA), col=c('black','red','blue'), 
    lty=c(NA,NA,1),lwd=c(NA,NA,2),bty='n',cex=0.75,title=" % p-values>0.05")
  
  #restore mar to default value 
  par(mar=c(5, 4, 4, 2) + 0.1)      
  dev.off()
  #bonferroni correction1
  unt_pval5<-round(length(merge[merge[['grp']]=='untransformed' & merge[['value']]>(0.05/649),][['value']])*100/
                     length(merge[merge[['grp']]=='untransformed',][['value']]),3) #13.885%
  unt_pval5*650*649*0.01*0.5 #29286.94
  #bonferroni2
  unt_pval5<-round(length(merge[merge[['grp']]=='untransformed' & merge[['value']]>(0.05/649*650*0.5),][['value']])*100/
                     length(merge[merge[['grp']]=='untransformed',][['value']]),3) #0.547%
  unt_pval5*650*649*0.01*0.5 #1153.76
}

