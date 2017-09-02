## try http:// if https:// URLs are not supported
#
# http://www.bioconductor.org/packages/2.12/data/experiment/html/HD2013SGI.html
#source("https://bioconductor.org/biocLite.R")
#biocLite("HD2013SGI")
#biocLite("LMGene") 

library(LMGene)
library(Biobase) 
library(tools)

# LWafula (25-06-2017) MRMRFS script practise
#library(MRMRFS)
library(hexbin)
library(ggplot2)
library(matlab)
library(plyr)
library(dplyr)
library(data.table)
library(flux)

#install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival") 
# source("http://bioconductor.org/biocLite.R") 
# biocLite(c("GO.db", "preprocessCore", "impute")) 

library(WGCNA) #https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/InstallationInstructions.html
library(knitr) #for kable
library("reshape2")


#to view documentation for the package
browseVignettes("HD2013SGI")
browseVignettes("LMGene")
browseVignettes("Biobase")

#analysis of Microarray data and glog-transformation
# Load the sample ExpressionSet class data in the package LMGene.
data(sample.eS)

#View the data structure of the sample data and the details of exprs and 
#phenoData slots in the data.
slotNames(sample.eS)
dim(exprs(sample.eS))
exprs(sample.eS)[1:3,]

#
phenoData(sample.eS)
slotNames(phenoData(sample.eS))

#LMGene User’s Guide
#Data generation. If you don’t have ExpressionSet class data, you need to make 
#some. LMGene provides a function that can generate an object of class 
#ExpressionSet, assuming that there are array data of matrix class and 
#experimental data of list class

data(sample.mat) 
dim(sample.mat)

data(vlist)
vlist

# Generate ExpressionSet class data using neweS function.
test.eS<-neweS(sample.mat, vlist) 
class(test.eS)

?neweS 
# CODES for access and some data ::  http://dmrocke.ucdavis.edu/software.html
# documentation ::  http://127.0.0.1:26830/library/LMGene/doc/LMGene.pdf

#GLOG parameter estimation
tranpar <- tranest(sample.eS)
tranpar

#controlling number of genes used in the parameter estimation
tranpar <- tranest(sample.eS, ngenes=100)
tranpar

#2. G-log transformation. Using the obtained two parameters, the g-log 
# transformed expression set can be calculated as follows.
trsample.eS <- transeS(sample.eS, tranpar$lambda, tranpar$alpha)
exprs(sample.eS)[1:3,1:8]

# 3. Tranest options: multiple alpha, lowessnorm, model
# Rather than using a single alpha for all samples, we can estimate a separate 
# alpha for each sample. This allows for differences in chips, in sample 
# concentration, or exposure conditions.
tranparmult <- tranest(sample.eS, mult=TRUE)
tranparmult

trsample.eS <- transeS (sample.eS, tranparmult$lambda, tranparmult$alpha)
exprs(trsample.eS)[1:3,1:8]


## TRY thz when everything is settled later on
load(gloglambda0_p1601xxPPS_means) #==>outfilelog0 #aggregated dataset for data1

#save a transposed version since the package expects rows to b features and 
# columns as the samples
gloglambda0_p1601xxPPS_meanscsv=
"E:/uhasselt/22yrsem/MThesis/data/gloglambda0_p1601xxPPS_meanscsv"

write.csv(outfilelog0,gloglambda0_p1601xxPPS_meanscsv,sep=',')

#hizo NAs and -Inf zitoe b4 ufike hapa, bure itakuwa noma
exprss=
as.matrix(read.table(gloglambda0_p1601xxPPS_meanscsv,sep=',',header=T),
,sep=',',header=T, row.names=1,as.is=T)

exprss=t(exprss)
minimalSet <- ExpressionSet(assayData=exprss)

#trying out the lambdas and alphas for glog dataset
slotNames(minimalSet)
dim(exprs(minimalSet))

#thz should be done in the dataPrep however!
tranpar <- tranest(minimalSet)
tranpar

tranparmult0 <- tranest(minimalSet, mult=TRUE)
tranparmult0



