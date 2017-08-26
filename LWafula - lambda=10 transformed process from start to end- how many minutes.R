

# LWafula 09062017

#this results are to be used to add to the untranformed plot for the Hotelling's results and see if there's a better sepration via the 
#distribution of the Hoteling's T-sqr statistic

# alpha estimates are available from the LWafula GLOG parameters estimation script in the MThesis/scripts folder
# the y-alpha dataset is also built in the LWafula GLOG parameters estimation script

# -here we are just combining the y-alpha datasets for every plate in a given datasetbatch
# -then tranform using proposal lambdas, investigate the time taken for a given lambda value to the end
# create a loop and see the results from this loop, we are saving the NFeatures selected, active calls

#lambda used here is lambda=10
# also did for lambda=25
# 
rm(list=ls())
tstart <- Sys.time() #start time

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

datawd<-'E:/uhasselt/22yrsem/MThesis/data'
dataglog<-'E:/uhasselt/22yrsem/MThesis/data/glog_estimates'

# add plate layout to know the controls
fn_in_wellinfo  <- 'p1601xxPPS_wellInfo_forExternalUse.csv'
#
fn_in_wellinfo  <- file.path(datawd, fn_in_wellinfo)

winfo <- read.csv(fn_in_wellinfo, as.is = TRUE) #14688x13



########## ``````````````````````````````````````````
plate67070<-'phaedra_singlecell_p160118PPS008__plate67070.Rdata'
plate67070<-file.path(datawd,plate67070)
load(plate67070)
plate67070<-dc

# 
plate67070 <- data.table(plate67070) 

#Identify the feature columns and remove troublesome features
cn0 <- names(plate67070)

#grepl returns TRUE if a string contains the pattern, otherwise FALSE
features0 <- cn0[grepl('^(Cell_|Cyto_|Nuc_|Ratio_|CellCount|CellConfluency)',
                       cn0, ig = TRUE)]

rm(list=c('plate67070','dc'))
# 
#rejoin the y-alpha datasets and run for a given lambda, save the resulting combined dataset for this y-alpha 

#lambda loop
#

  # batch 1
  
  #batch1              # for lambda =0.1
  file.names <- dir(dataglog, pattern ="glogparams_batch1_")
  gloglambda10<-''
  
  for(i in 1:length(file.names)){
    fileb1<-file.path(dataglog,file.names[i])
    load(fileb1)
    # lambda adjusts
    fileb1[,features0]<-log(fileb1[,features0]+sqrt(fileb1[,features0]^2 +10)); fileb1<-data.table(fileb1)
    featuremeansbywell=fileb1[, c(lapply(.SD, function(x) mean(x[!is.infinite(x)], na.rm = TRUE)), .(CellCount_AllRetained = .N)),
                              by= .(WELL_ID, ROW_NR, COL_NR, PLATE_ID), .SDcols = features0]
    gloglambda10 <- rbind(gloglambda10, featuremeansbywell,fill=T)
    rm(list=c('fileb1','featuremeansbywell'))
  }
  gloglambda10=gloglambda10[-c(1),-1]
  
  #batch2             # for lambda =0.1
  file.names <- dir(dataglog, pattern ="glogparams_batch2_")
  for(i in 1:length(file.names)){
    fileb1<-file.path(dataglog,file.names[i])
    load(fileb1)
    # lambda adjusts
    fileb1[,features0]<-log(fileb1[,features0]+sqrt(fileb1[,features0]^2 +10)); fileb1<-data.table(fileb1)
    featuremeansbywell=fileb1[, c(lapply(.SD, function(x) mean(x[!is.infinite(x)], na.rm = TRUE)), .(CellCount_AllRetained = .N)),
                              by= .(WELL_ID, ROW_NR, COL_NR, PLATE_ID), .SDcols = features0]
    gloglambda10 <- rbind(gloglambda10, featuremeansbywell,fill=T)
    rm(list=c('fileb1','featuremeansbywell'))
  }
  
  #batch3             # for lambda =0.1
  file.names <- dir(dataglog, pattern ="glogparams_batch3_")
  for(i in 1:length(file.names)){
    fileb1<-file.path(dataglog,file.names[i])
    load(fileb1)
    # lambda adjusts
    fileb1[,features0]<-log(fileb1[,features0]+sqrt(fileb1[,features0]^2 +10)); fileb1<-data.table(fileb1)
    featuremeansbywell=fileb1[, c(lapply(.SD, function(x) mean(x[!is.infinite(x)], na.rm = TRUE)), .(CellCount_AllRetained = .N)),
                              by= .(WELL_ID, ROW_NR, COL_NR, PLATE_ID), .SDcols = features0]
    gloglambda10 <- rbind(gloglambda10, featuremeansbywell,fill=T)
    rm(list=c('fileb1','featuremeansbywell'))
  }
  
  ## then do the data preprocessing as usual
  # data dimension before aggregating for dataset 1
  
  
  #################### ```````````````````````````````` ###############################################################################
  ##          `````````````````````````DATA PREPROCESSING `````````````````````````````````````````````````````````````````````````##`
  
  # normalization and convenience functions
  source('E:/uhasselt/22yrsem/MThesis/scripts/R_code/stools.R')
  
  #WD
  genwd<-'E:/uhasselt/22yrsem/MThesis'
  datadir<-'E:/uhasselt/22yrsem/MThesis/data'
  dataview<-'E:/uhasselt/22yrsem/MThesis/report/dataView'
  # # combined data for labda=j 
  str(gloglambda10, list.len=ncol(gloglambda10))
  
  ## ``` Add plate layout (compound, concentration, .) ````
  
  glogdatp10 <- 
    merge(gloglambda10, winfo, by = c('WELL_ID', 'PLATE_ID','COL_NR','ROW_NR'), 
          suffixes = c('', '.from_well_info'))
  stopifnot(!any(is.na(glogdatp10$Treatment)))
  
  # count number of replicates for each treatment over wells
  glogdatp10[, n := .N, by = Treatment]
  table(glogdatp10[['Treatment']],glogdatp10[['n']]); 
  dim(table(glogdatp10[['Treatment']],glogdatp10[['n']])) #1253 compounds
  
  #PLATES (both barcode and plateID are plate IDs) #plates have 248/260/268/272/288/296 wells
  sort(table(glogdatp10$BARCODE))
  
  #number of replicates for the samples
  # View(glogdatp10[,c('WELLTYPE_CODE','n','Treatment')])
  n <- glogdatp10[WELLTYPE_CODE == 'SAMPLE', n[1], by = Treatment]$V1
  
  #hist(n, main = 'Number of replicates', xlab = 'Number of Replicates', 
  #ylab = 'Number of Treatments', breaks = 1 : max(n))
  
  ## code below assumes data.frame (not data.table)
  glogdatp10<- as.data.frame(glogdatp10) 
  
  #Identify the feature columns and remove troublesome features
  cn <- names(glogdatp10)
  
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
  st <- data.frame.OverviewColumnStats.2(glogdatp10)
  st <- arrange(st, -n_NA_Inf_0)
  browse_table(st, 'stats input_all columns')
  kable(st[st$n_NA_Inf_0 > 0, ], caption = 'Column stats (all columns)')
  
  #features only
  st <- data.frame.OverviewColumnStats.2(glogdatp10[, features])
  st <- arrange(st, -n_NA_Inf_0)
  browse_table(st, 'stats input_features only')
  kable(st[st$n_NA_Inf_0 > 0, ], caption = 'Column stats (features only)')
  
  # at normalization, we've 461 features, 14688 wellIDs
  # Normalized
  cellcountfeatures_raw <- colnames(glogdatp10)[grepl('CellCount_AllRetained', colnames(glogdatp10))]
  stopifnot(length(cellcountfeatures_raw) > 0)
  
  # NOTE; check later that cellcountfeatures is sorted
  
  glogdatp10<- normalize.percentControls(
    glogdatp10, 
    groupBy             = 'BARCODE',
    WELLTYPE_CODES          = c('LC'),                          # WELLTYPE_CODE(S) of the control wells
    featureColumns      = cellcountfeatures_raw,      # columns to normalize, NULL: all numeric columns not in 'colnames.ignore'
    agg.function        = median,
    suffix              = '.prctCntr',
    colnames.ignore     = NULL
  )
  
  # Normalize all features: z-score relative to DMSO controls
  glogdatp10<- normalize.zscoreControls_plateMean_globalPooledStdev(
    glogdatp10,
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
  features_norm <- attr(glogdatp10, 'features.normalized')
  
  # drop raw features from data frame, except cell count feature
  features_raw_drop <- features[!grepl('CellCount_AllRetained', features)]
  glogdatp10<- glogdatp10[, setdiff(names(glogdatp10), features_raw_drop)]
  
  cat(paste0(sort(features_norm), collapse = '\n'))
  
  # Column statistics (after normalization)
  #check that there are no NAs
  st <- data.frame.OverviewColumnStats.2(glogdatp10[, features_norm])
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
  valueControl <- unique(glogdatp10[glogdatp10$WELLTYPE_CODE == 'LC', 
                                  columnReplicateID])
  stopifnot(length(valueControl) == 1)
  
  # ----------
  #  run mRMR
  # ----------
  FSres10 <- run_mRMR(
    dataset = glogdatp10,
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
  # character vector features_norm_sel #19 features selected
  features_norm_sel <- FSres10$selectedFeatures[1:FSres10$optN] # select the 'optimal' number of features
  lambdafeatactv_alphacorrected[j,1]<-10
  lambdafeatactv_alphacorrected[j,2]<-length(features_norm_sel)
  
  kable(data.frame(Rank = 1 : length(features_norm_sel), Feature = features_norm_sel))
  
  stopifnot(!any(is.na(glogdatp10$Treatment)))
  
  ## ```` In/active calling ````
  ## Treatments that don't show a significant difference to DMSO controls are flagged as inactive. 
  ## This information can be used in downstream analysis to work with only the active treatments where 
  ## appropriate.
  ## % active treatments:  56.42458 meaning that 56.42458% of the CPDs have a significant different to DMSO
  
  source('E:/uhasselt/22yrsem/MThesis/scripts/R_code/active_calling.R')
  glogdatp10<- active_calling(glogdatp10, features_norm_sel, active_distance_percentile = 0.95)
  
  stopifnot(!any(is.na(glogdatp10$Treatment))) 

  #active treatments. selected features #
  datpfeatact10<-data.table(glogdatp10)
  datpfeatact10<-datpfeatact10[Active_FractionOfRepl>0.5,.SD, .SDcols=c(features_norm_sel, 'Treatment','CMPD_ID','CONCENTRATION')]
  
  fn_out10          <- 'p1601xxPPS_means_preprocessedfeat_actlog10.Rdata'
  fn_out10          <- file.path(datadir, fn_out10)
  save(datpfeatact10, file = fn_out10)
  cat('written to ',fn_out10,'\n')
