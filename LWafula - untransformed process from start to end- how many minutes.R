
#untransformed took 13 minutes
# for the untransormed and transformed dataset, how long does it take to run the code to the end?
# LWafula: 06/06/2017
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

datawd<-'E:/uhasselt/22yrsem/MThesis/data'
plate67070<-'phaedra_singlecell_p160118PPS008__plate67070.Rdata'

plate67070<-file.path(datawd,plate67070)
load(plate67070)
plate67070<-dc

#remove problematic features 
plate67070 <- data.table(plate67070) 

#Identify the feature columns and remove troublesome features
cn0 <- names(plate67070)

#grepl returns TRUE if a string contains the pattern, otherwise FALSE
features0 <- cn0[grepl('^(Cell_|Cyto_|Nuc_|Ratio_|CellCount|CellConfluency)',
                       cn0, ig = TRUE)]

rm(list=c('plate67070','dc'))

#plate layout data
fn_in_wellinfo  <- 'p1601xxPPS_wellInfo_forExternalUse.csv'
fn_in_wellinfo  <- file.path(datawd, fn_in_wellinfo)

winfo <- read.csv(fn_in_wellinfo, as.is = TRUE) #14688x13


# UNTRANSFORMED
  # batch 1
  file.names <- dir(datawd, pattern ="phaedra_singlecell_p160118PPS008__plate")
  outfilelog1<-''
  
  for(i in 1:length(file.names)){
    file<-file.path(datawd,file.names[i])
    load(file)
    file <-dc
    featuremeansbywell=file[, c(lapply(.SD, function(x) mean(x[!is.infinite(x)], na.rm = TRUE)), .(CellCount_AllRetained = .N)),
                            by= .(WELL_ID, ROW_NR, COL_NR, PLATE_ID), .SDcols = features0]
    outfilelog1 <- rbind(outfilelog1, featuremeansbywell,fill=T)
    rm(list=c('file','dc','featuremeansbywell'))
  }
  outfilelog1=outfilelog1[complete.cases(outfilelog1$PLATE_ID),-c('x')]
  
  #batch 2
  file.names <- dir(datawd, pattern ="phaedra_singlecell_p160125PPS009__plate")
  for(i in 1:length(file.names)){
    file<-file.path(datawd,file.names[i])
    load(file)
    file <-dc
    featuremeansbywell=file[, c(lapply(.SD, function(x) mean(x[!is.infinite(x)], na.rm = TRUE)), .(CellCount_AllRetained = .N)),
                            by= .(WELL_ID, ROW_NR, COL_NR, PLATE_ID), .SDcols = features0]
    outfilelog1 <- rbind(outfilelog1, featuremeansbywell,fill=T)
    rm(list=c('file','featuremeansbywell','dc'))
  }
  
  #batch 3
  file.names <- dir(datawd, pattern ="phaedra_singlecell_p160201PPS010__plate")
  for(i in 1:length(file.names)){
    file<-file.path(datawd,file.names[i])
    load(file)
    file <-dc
    featuremeansbywell=file[, c(lapply(.SD, function(x) mean(x[!is.infinite(x)], na.rm = TRUE)), .(CellCount_AllRetained = .N)),
                            by= .(WELL_ID, ROW_NR, COL_NR, PLATE_ID), .SDcols = features0]
    outfilelog1 <- rbind(outfilelog1, featuremeansbywell,fill=T)
    rm(list=c('file','featuremeansbywell','dc'))
  }
  
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
  str(outfilelog1, list.len=ncol(outfilelog1))
  
  ## ``` Add plate layout (compound, concentration, .) ````
  
  datp <- 
    merge(outfilelog1, winfo, by = c('WELL_ID', 'PLATE_ID','COL_NR','ROW_NR'), 
          suffixes = c('', '.from_well_info'))
  stopifnot(!any(is.na(datp$Treatment)))
  
  # count number of replicates for each treatment over wells
  datp[, n := .N, by = Treatment]
  table(datp[['Treatment']],datp[['n']]); 
  dim(table(datp[['Treatment']],datp[['n']])) #1253 compounds
  
  #PLATES (both barcode and plateID are plate IDs) #plates have 248/260/268/272/288/296 wells
  sort(table(datp$BARCODE))
  
  #number of replicates for the samples
  # View(datp[,c('WELLTYPE_CODE','n','Treatment')])
  n <- datp[WELLTYPE_CODE == 'SAMPLE', n[1], by = Treatment]$V1
  
  #hist(n, main = 'Number of replicates', xlab = 'Number of Replicates', 
  #ylab = 'Number of Treatments', breaks = 1 : max(n))
  
  ## code below assumes data.frame (not data.table)
  datp<- as.data.frame(datp) 
  
  #Identify the feature columns and remove troublesome features
  cn <- names(datp)
  
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
  st <- data.frame.OverviewColumnStats.2(datp)
  st <- arrange(st, -n_NA_Inf_0)
  browse_table(st, 'stats input_all columns')
  kable(st[st$n_NA_Inf_0 > 0, ], caption = 'Column stats (all columns)')
  
  #features only
  st <- data.frame.OverviewColumnStats.2(datp[, features])
  st <- arrange(st, -n_NA_Inf_0)
  browse_table(st, 'stats input_features only')
  kable(st[st$n_NA_Inf_0 > 0, ], caption = 'Column stats (features only)')
  
  # at normalization, we've 461 features, 14688 wellIDs
  # Normalized
  cellcountfeatures_raw <- colnames(datp)[grepl('CellCount_AllRetained', colnames(datp))]
  stopifnot(length(cellcountfeatures_raw) > 0)
  
  # NOTE; check later that cellcountfeatures is sorted
  
  datp<- normalize.percentControls(
    datp, 
    groupBy             = 'BARCODE',
    WELLTYPE_CODES          = c('LC'),                          # WELLTYPE_CODE(S) of the control wells
    featureColumns      = cellcountfeatures_raw,      # columns to normalize, NULL: all numeric columns not in 'colnames.ignore'
    agg.function        = median,
    suffix              = '.prctCntr',
    colnames.ignore     = NULL
  )
  
  # Normalize all features: z-score relative to DMSO controls
  datp<- normalize.zscoreControls_plateMean_globalPooledStdev(
    datp,
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
  features_norm <- attr(datp, 'features.normalized')
  
  # drop raw features from data frame, except cell count feature
  features_raw_drop <- features[!grepl('CellCount_AllRetained', features)]
  datp<- datp[, setdiff(names(datp), features_raw_drop)]
  
  cat(paste0(sort(features_norm), collapse = '\n'))
  
  # Column statistics (after normalization)
  #check that there are no NAs
  st <- data.frame.OverviewColumnStats.2(datp[, features_norm])
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
  valueControl <- unique(datp[datp$WELLTYPE_CODE == 'LC', 
                                        columnReplicateID])
  stopifnot(length(valueControl) == 1)
  
  # ----------
  #  run mRMR
  # ----------
  FSres <- run_mRMR(
    dataset = datp,
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
  features_norm_sel <- FSres$selectedFeatures[1:FSres$optN] # select the 'optimal' number of features
  
  kable(data.frame(Rank = 1 : length(features_norm_sel), Feature = features_norm_sel))
  
  stopifnot(!any(is.na(datp$Treatment)))
  
  ## ```` In/active calling ````
  ## Treatments that don't show a significant difference to DMSO controls are flagged as inactive. 
  ## This information can be used in downstream analysis to work with only the active treatments where 
  ## appropriate.
  ## % active treatments:  56.42458 meaning that 56.42458% of the CPDs have a significant different to DMSO
  
  source('E:/uhasselt/22yrsem/MThesis/scripts/R_code/active_calling.R')
  datp<- active_calling(datp, features_norm_sel, active_distance_percentile = 0.95)
  
 dim(table(datp[datp$Active_FractionOfRepl>0.5,]$Treatment))*100/1253
  
  stopifnot(!any(is.na(datp$Treatment))) 
  
  # transormed one
  
  




