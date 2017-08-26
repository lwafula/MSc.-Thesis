# DATASET 2
# meeting held on 06062017 deliberations

# alpha estimates are available from the LWafula GLOG parameters estimation script in the MThesis/DATASET2/scripts folder
# the y-alpha dataset is also built in the LWafula GLOG parameters estimation script

# -here we are just combining the y-alpha datasets for every plate in a given datasetbatch
# -then tranform using proposal lambdas, investigate the time taken for a given lambda value to the end
# create a loop and see the results from this loop, we are saving the NFeatures selected, active calls
# save datasets for each lambda

rm(list=ls())
tstart <- Sys.time() #start time

setwd('E:/uhasselt/22yrsem/MThesis/data')
library(hexbin)
library(ggplot2)
library(matlab)
library(plyr)
library(dplyr)
library(data.table)
library(flux)
library(WGCNA)
library(knitr) #for kable

datawd<-'E:/uhasselt/22yrsem/MThesis/data'
dataglog<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog_estimates'
dataglog_data<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda'

# add plate layout to know the controls
fn_in_wellinfo  <- 'p1602xxPPS_wellInfo_forExternalUse.csv'
#
fn_in_wellinfo  <- file.path(datawd, fn_in_wellinfo)

winfo <- read.csv(fn_in_wellinfo, as.is = TRUE) #14688x13



########## ``````````````````````````````````````````
plate67413<-'phaedra_singlecell_p160208PPS011__plate67413.Rdata'
plate67413<-file.path(datawd,plate67413)
load(plate67413)
plate67413<-dc

# 
plate67413 <- data.table(plate67413) 

#Identify the feature columns and remove troublesome features
cn0 <- names(plate67413)

#grepl returns TRUE if a string contains the pattern, otherwise FALSE
features0 <- cn0[grepl('^(Cell_|Cyto_|Nuc_|Ratio_|CellCount|CellConfluency)',
                       cn0, ig = TRUE)]

rm(list=c('plate67413','dc'))
# 
#rejoin the y-alpha datasets and run for a given lambda, save the resulting combined dataset for this y-alpha 


# @@@@@@@@@@@@@@@@@@@ RUN THZ LOOP WHEN U ABSOLUTELY WANT TO @@@@@@@@@@@@@@@@@@@@@@@@@@@ #
#lambda loop j=50 is lambda=25,j=21 is lambda=10.5
#
lambda=seq(0.1,15,by=0.1)
lambdafeatactv_alphacorrected=array(NA,dim=c(length(lambda),3))

for(j in 1:length(lambda)){
  # batch 1
  
  #batch1              # for lambda =0.1
  file.names <- dir(dataglog, pattern ="glogparams_batch1_")
  gloglambda<-''
  
  for(i in 1:length(file.names)){
    fileb1<-file.path(dataglog,file.names[i])
    load(fileb1)
    # lambda adjusts
    fileb1[,features0]<-log(fileb1[,features0]+sqrt(fileb1[,features0]^2 +lambda[j])); fileb1<-data.table(fileb1)
    featuremeansbywell=fileb1[, c(lapply(.SD, function(x) mean(x[!is.infinite(x)], na.rm = TRUE)), .(CellCount_AllRetained = .N)),
                              by= .(WELL_ID, ROW_NR, COL_NR, PLATE_ID), .SDcols = features0]
    gloglambda <- rbind(gloglambda, featuremeansbywell,fill=T)
    rm(list=c('fileb1','featuremeansbywell'))
  }
  gloglambda=gloglambda[-c(1),-1]
  
  #batch2             # for lambda =0.1
  file.names <- dir(dataglog, pattern ="glogparams_batch2_")
  for(i in 1:length(file.names)){
    fileb1<-file.path(dataglog,file.names[i])
    load(fileb1)
    # lambda adjusts
    fileb1[,features0]<-log(fileb1[,features0]+sqrt(fileb1[,features0]^2 +lambda[j])); fileb1<-data.table(fileb1)
    featuremeansbywell=fileb1[, c(lapply(.SD, function(x) mean(x[!is.infinite(x)], na.rm = TRUE)), .(CellCount_AllRetained = .N)),
                              by= .(WELL_ID, ROW_NR, COL_NR, PLATE_ID), .SDcols = features0]
    gloglambda <- rbind(gloglambda, featuremeansbywell,fill=T)
    rm(list=c('fileb1','featuremeansbywell'))
  }
  
  #batch3             # for lambda =0.1
  file.names <- dir(dataglog, pattern ="glogparams_batch3_")
  for(i in 1:length(file.names)){
    fileb1<-file.path(dataglog,file.names[i])
    load(fileb1)
    # lambda adjusts
    fileb1[,features0]<-log(fileb1[,features0]+sqrt(fileb1[,features0]^2 +lambda[j])); fileb1<-data.table(fileb1)
    featuremeansbywell=fileb1[, c(lapply(.SD, function(x) mean(x[!is.infinite(x)], na.rm = TRUE)), .(CellCount_AllRetained = .N)),
                              by= .(WELL_ID, ROW_NR, COL_NR, PLATE_ID), .SDcols = features0]
    gloglambda <- rbind(gloglambda, featuremeansbywell,fill=T)
    rm(list=c('fileb1','featuremeansbywell'))
  }
  gloglambda.temp<-gloglambda
  # column identifying the lambdavalue
  
  
  ## then do the data preprocessing as usual
  # data dimension before aggregating for dataset 1
  
  
  #################### ```````````````````````````````` ###############################################################################
  ##          `````````````````````````DATA PREPROCESSING `````````````````````````````````````````````````````````````````````````##`
  
  # normalization and convenience functions
  source('E:/uhasselt/22yrsem/MThesis/scripts/R_code/stools.R')
  
  #WD
  genwd<-'E:/uhasselt/22yrsem/MThesis/DATASET2'
  datadir<-'E:/uhasselt/22yrsem/MThesis/data'
  dataview<-'E:/uhasselt/22yrsem/MThesis/DATASET2/report/dataView'
  # # combined data for labda=j 
  str(gloglambda, list.len=ncol(gloglambda))
  
  ## ``` Add plate layout (compound, concentration, .) ````
  
  glogdatp <- 
    merge(gloglambda, winfo, by = c('WELL_ID', 'PLATE_ID','COL_NR','ROW_NR'), 
          suffixes = c('', '.from_well_info'))
  stopifnot(!any(is.na(glogdatp$Treatment)))
  
  # count number of replicates for each treatment over wells
  glogdatp[, n := .N, by = Treatment]
  table(glogdatp[['Treatment']],glogdatp[['n']]); 
  dim(table(glogdatp[['Treatment']],glogdatp[['n']])) #1253 compounds
  
  #PLATES (both barcode and plateID are plate IDs) #plates have 248/260/268/272/288/296 wells
  sort(table(glogdatp$BARCODE))
  
  #number of replicates for the samples
  # View(glogdatp[,c('WELLTYPE_CODE','n','Treatment')])
  n <- glogdatp[WELLTYPE_CODE == 'SAMPLE', n[1], by = Treatment]$V1
  
  #hist(n, main = 'Number of replicates', xlab = 'Number of Replicates', 
  #ylab = 'Number of Treatments', breaks = 1 : max(n))
  
  ## code below assumes data.frame (not data.table)
  glogdatp<- as.data.frame(glogdatp) 
  
  #Identify the feature columns and remove troublesome features
  cn <- names(glogdatp)
  
  #grepl returns TRUE if a string contains the pattern, otherwise FALSE
  features <- cn[grepl('^(Cell_|Cyto_|Nuc_|Ratio_|CellCount|CellConfluency)',
                       cn, ig = TRUE)]
  
  #notice the non-features
  setdiff(cn,features)
  
  features <- features[!grepl('Nuc_.*_SER.*_4', features)]      # remove texture features at scale 4 for nuclei (small region ==> NaN)
  features <- features[!grepl('NearestNeighbourDist', features)]              # too many zeros
  features <- features[!grepl('Cyto_NucCytoBFP_SERSpot_4.median', features)]  # too many zeros
  features <- features[!grepl('Cyto_NucCytoBFP_SERSpot_2.median', features)]  # too many zeros
  features <- features[!grepl('CellCount_AllDetected', features)]             # CellCount_AllDetected is redundant, CellCount_AllRetained contains the actual number of nuclei/cells in the well (after removing those touching the image borders)
  
  # features <- features[sample(1:length(features), size = 100)] # speed-up for debugging!!!
  
  cat('Number of features: ', length(features), '\n') #461 features from 463
  
  # ```column statistics before normalization```
  st <- data.frame.OverviewColumnStats.2(glogdatp)
  st <- arrange(st, -n_NA_Inf_0)
  browse_table(st, 'stats input_all columns')
  kable(st[st$n_NA_Inf_0 > 0, ], caption = 'Column stats (all columns)')
  
  #features only
  st <- data.frame.OverviewColumnStats.2(glogdatp[, features])
  st <- arrange(st, -n_NA_Inf_0)
  browse_table(st, 'stats input_features only')
  kable(st[st$n_NA_Inf_0 > 0, ], caption = 'Column stats (features only)')
  
  # at normalization, we've 461 features, 14688 wellIDs
  # Normalized
  cellcountfeatures_raw <- colnames(glogdatp)[grepl('CellCount_AllRetained', colnames(glogdatp))]
  stopifnot(length(cellcountfeatures_raw) > 0)
  
  # NOTE; check later that cellcountfeatures is sorted
  
  glogdatp<- normalize.percentControls(
    glogdatp, 
    groupBy             = 'BARCODE',
    WELLTYPE_CODES          = c('LC'),                          # WELLTYPE_CODE(S) of the control wells
    featureColumns      = cellcountfeatures_raw,      # columns to normalize, NULL: all numeric columns not in 'colnames.ignore'
    agg.function        = median,
    suffix              = '.prctCntr',
    colnames.ignore     = NULL
  )
  
  # Normalize all features: z-score relative to DMSO controls
  glogdatp<- normalize.zscoreControls_plateMean_globalPooledStdev(
    glogdatp,
    groupBy                       = c('expID'),                   # these groups will be handled completely independently
    groupBy.mean            = c('BARCODE'),               # grouping for means (or medians)
    groupBy.pooled_sd           = NULL,                           # grouping for s.d. (or mad), pool the s.d. (or mad) into one value, NULL: across all plates
    WELLTYPE_CODES                = c('LC'),                  # WELLTYPE_CODE(S) of the control wells
    featureColumns          = features,                   # columns to normalize, NULL: all numeric columns not in 'colnames.ignore'
    suffix                  = '.zGSLC',
    centerBeforeComputingSD = TRUE,                             # TRUE: subtract within each "groupBy.mean"-group the mean (or median), then compute s.d. // FALSE: compute s.d. on non-centered data
    robust                  = FALSE,                  # FALSE: mean / sd, TRUE: median / mad
    colnames.ignore           = NULL
  )
  features_norm <- attr(glogdatp, 'features.normalized')
  
  # drop raw features from data frame, except cell count feature
  features_raw_drop <- features[!grepl('CellCount_AllRetained', features)]
  glogdatp<- glogdatp[, setdiff(names(glogdatp), features_raw_drop)]
  
  cat(paste0(sort(features_norm), collapse = '\n'))
  
  # Column statistics (after normalization)
  #check that there are no NAs
  st <- data.frame.OverviewColumnStats.2(glogdatp[, features_norm])
  st <- arrange(st, -n_NA_Inf_0)
  browse_table(st, 'stats after normalization_features only')
  kable(st[st$n_NA_Inf_0 > 0, ], caption = 'Column stats after normalization (features only)')
  
  # ````````````` Feature Selection ```````````````````````
  #loading package MRMRFs
  library(hcs)
  library(MRMRFS)
  
  source('E:/uhasselt/22yrsem/MThesis/scripts/R_code/featureSelection_mRMR.R')
  
  
  columnReplicateID   <- "Treatment"  # the column that contains the replicate information; all rows with the same value are considered replicates
  cellcount_ft            <- 'CellCount_AllRetained'
  min_cellcount       <- 100
  nFeatures               <- 75   # number of features to select (initially, the optimal number which is usually less than nFeatures will be determined below)
  N_REPLICATES            <- 0    # > 0: randomly sample n replicates per treatment (to reduce influence of compounds with large number of replicates)
  AUCeval.m <- 10
  AUCeval.nPairsReplicates <- 0.2 
  AUCeval.nPairsNonReplicates <- 0.2
  AUCeval.nPairsReplicatesMax <- 10000
  AUCeval.nPairsNonReplicatesMax <- 10000
  AUCeval.distMethod <- "pearson"
  AUCeval.nThresholdsAUC <- 100
  
  DIR_OUT <- 'C:/TEMP'
  
  # -------------------------------
  #  specify the control value
  # -------------------------------
  valueControl <- unique(glogdatp[glogdatp$WELLTYPE_CODE == 'LC', 
                                  columnReplicateID])
  stopifnot(length(valueControl) == 1)
  
  # ----------
  #  run mRMR
  # ----------
  FSres <- run_mRMR(
    dataset = glogdatp,
    features = features_norm,    # all features to select from
    nFeatures = nFeatures,            # number of features to be selected
    min_cellcount = min_cellcount,        # select features based on wells with at least this many cells
    cellcount_ft = cellcount_ft, 
    N_REPLICATES = N_REPLICATES, # > 0: randomly sample n replicates per treatment (to reduce influence of compounds with large number of replicates)
    columnReplicateID = columnReplicateID,
    valueControl = valueControl,
    AUCeval.m = AUCeval.m,
    AUCeval.nPairsReplicates = AUCeval.nPairsReplicates, 
    AUCeval.nPairsNonReplicates = AUCeval.nPairsNonReplicates,
    AUCeval.nPairsReplicatesMax = AUCeval.nPairsReplicatesMax,
    AUCeval.nPairsNonReplicatesMax = AUCeval.nPairsNonReplicatesMax,
    AUCeval.distMethod = AUCeval.distMethod,
    AUCeval.nThresholdsAUC = AUCeval.nThresholdsAUC,
    DIR_OUT = DIR_OUT
  )
  
  # The final feature set that makes up the pheno-signature of each sample is saved in the 
  # character vector features_norm_sel #38 features selected
  features_norm_sel <- FSres$selectedFeatures[1:FSres$optN] # select the 'optimal' number of features
  lambdafeatactv_alphacorrected[j,1]<-lambda[j]
  lambdafeatactv_alphacorrected[j,2]<-length(features_norm_sel)
  
  kable(data.frame(Rank = 1 : length(features_norm_sel), Feature = features_norm_sel))
  
  stopifnot(!any(is.na(glogdatp$Treatment)))
  
  ## ```` In/active calling ````
  ## Treatments that don't show a significant difference to DMSO controls are flagged as inactive. 
  ## This information can be used in downstream analysis to work with only the active treatments where 
  ## appropriate.
  ## % active treatments:  52.67358 meaning that 52.67358% of the CPDs have a significant different to DMSO
  
  source('E:/uhasselt/22yrsem/MThesis/scripts/R_code/active_calling.R')
  glogdatp<- active_calling(glogdatp, features_norm_sel, active_distance_percentile = 0.95)
  
  lambdafeatactv_alphacorrected[j,3]<-dim(table(glogdatp[glogdatp$Active_FractionOfRepl>=0.5,]$Treatment))*100/1253
  
  
  stopifnot(!any(is.na(glogdatp$Treatment)))
  
  # Save preprocessed data with selected features to file
  attr(glogdatp, 'features')           <- features_norm 
  attr(glogdatp, 'features_mRMR_all')  <- FSres$selectedFeatures
  attr(glogdatp, 'features_mRMR_optN') <- features_norm_sel
  
  #save 
  fn_out<- paste0('p1602xxPPS_means_preprocessedfeat_act_',lambda[j],'.RData')
  fn_out<-file.path(dataglog_data,fn_out)
  save(glogdatp, file = fn_out)
  #cat('written to ',fn_out0,'\n')
  print (lambda[j])
  ### imeisha hii #############
}



#lambdafeatactv_alphacorrected_temphold<-lambdafeatactv_alphacorrected
tend <- Sys.time() #end time

#time taken
c(tstart,tend,as.Date(tend,origin="2017-06-07")-as.Date(tstart,origin='2017-06-07'))

# @@@@@@@@@@@@@@@@@@@ RUN THZ LOOP WHEN U ABSOLUTELY WANT TO @@@@@@@@@@@@@@@@@@@@@@@@@@@ # 

# # #plots
# lambdavsoptfeature_alpha<-'D2_lambdavsoptfeature_alpha.pdf'
# lambdavsoptfeature_alpha<-file.path(dataview,lambdavsoptfeature_alpha)
# pdf(lambdavsoptfeature_alpha)
# plot(lambdafeatactv_alphacorrected[,1],lambdafeatactv_alphacorrected[,2],pch=19,type='l',xlab=expression(lambda),
#      ylab='Optimal features',ylim=c(25,45))
# abline(h=34, lty=2,lwd=1.5)
# legend('topright',legend=c(expression(paste('transformed by',' ',lambda,' ','=x-value')),'untransformed'),lty=c(1,2),bty='n')
# dev.off()
# # # 
# lambdavsactcall_alpha<-'D2_lambdavsactcall_alpha.pdf'
# lambdavsactcall_alpha<-file.path(dataview,lambdavsactcall_alpha)
# pdf(lambdavsactcall_alpha)
# plot(lambdafeatactv_alphacorrected[,1],lambdafeatactv_alphacorrected[,3],pch=19,type='l',xlab=expression(lambda),
#      ylab='active calling %',ylim=c(50,55))
# abline(h=51.23703, lty=2,lwd=1.5)
# legend('topright',legend=c(expression(paste('transformed by',' ',lambda,' ','=x-value')),'untransformed'),lty=c(1,2),bty='n')
# dev.off()
# 
# #save
lambdafeatactv_alphacorrected0<-'E:/uhasselt/22yrsem/MThesis/DATASET2/data/glog transformed datasets by lambda/lambdafeatactv_alphacorrected0p1to15.csv'
write.csv(lambdafeatactv_alphacorrected,lambdafeatactv_alphacorrected0)






