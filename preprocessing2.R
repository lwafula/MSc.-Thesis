# remove trouble-some features, normalize and select features
# same as preprocessing.R, but starting from the single-cell aggregated data (well info needs to be added) instead of the Phaedra well-level data

rm(list = ls())
setwd('C:/sjaensch/dev/R')

source('PhenoPrint/MasterThesis2017/stools.R') # normalization and convenience functions

datadir <- 'C:/sjaensch/docu/HCS projects/2014-09 ITS Proposal -- PhenoPrinting/MasterThesis 2017/data'


fn_in_means     <- 'p1602xxPPS_means.Rdata'
fn_in_wellinfo  <- 'p1602xxPPS_wellInfo_forExternalUse.csv'
fn_out          <- 'p1602xxPPS_means_preprocessed.Rdata'

fn_in_means     <- 'p1601xxPPS_means.Rdata'
fn_in_wellinfo  <- 'p1601xxPPS_wellInfo_forExternalUse.csv'
fn_out          <- 'p1601xxPPS_means_preprocessed.Rdata'

fn_in_means     <- file.path(datadir, fn_in_means)
fn_in_wellinfo  <- file.path(datadir, fn_in_wellinfo)
fn_out          <- file.path(datadir, fn_out)

load(fn_in_means) # ==> dw

#:::::::::::::::::::::::::::::
#
#  add well info
#
#:::::::::::::::::::::::::::::
# Note: The plate layout has been extracted from data from the PhenoPrint workflow. 
#       It contains only valid wells, in particular fluorescent compounds are already removed. 
#       When joining the well data with the plate layout, the wells not present in the plate layout 
#       will therefore automatically be dropped.

winfo <- read.csv(fn_in_wellinfo, as.is = TRUE)
datp <- merge(dw, winfo, by = c('WELL_ID', 'PLATE_ID', 'ROW_NR', 'COL_NR'), suffixes = c('', '.from_well_info')) # note: merging by WELL_ID alone would actually be sufficient as WELL_ID is a globally unique identifier
stopifnot(!any(is.na(datp$Treatment)))

#:::::::::::::::::::::::::::::
#
#  data.table --> data.frame
#
#:::::::::::::::::::::::::::::
datp <- as.data.frame(datp) # code below assumes data.frame (not data.table)

#:::::::::::::::::::::::::::::
#
#  feature columns
#
#:::::::::::::::::::::::::::::
cn <- names(datp)
features <- cn[grepl('^(Cell_|Cyto_|Nuc_|Ratio_|CellCount|CellConfluency)', cn, ig = TRUE)]
features <- features[!grepl('Nuc_.*_SER.*_4', features)]                    # remove texture features at scale 4 for nuclei (small region ==> NaN)
features <- features[!grepl('NearestNeighbourDist', features)]              # too many zeros
features <- features[!grepl('Cyto_NucCytoBFP_SERSpot_4.median', features)]  # too many zeros
features <- features[!grepl('Cyto_NucCytoBFP_SERSpot_2.median', features)]  # too many zeros
features <- features[!grepl('CellCount_AllDetected', features)]             # CellCount_AllRetained is sufficient, containing the kept number of nuclei/cells

cat(paste0(sort(features), collapse = '\n'))

#:::::::::::::::::::::::::::::
#
#  stats (zeros? NA?)
#
#:::::::::::::::::::::::::::::
browse_table(data.frame.OverviewColumnStats.2(datp), 'stats input')
browse_table(data.frame.OverviewColumnStats.2(datp[, features]), 'stats input _ features only')

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#  normalize cell count: % of DMSO controls
#
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cellcountfeatures_raw <- colnames(datp)[grepl('CellCount_AllRetained', colnames(datp))]

stopifnot(length(cellcountfeatures_raw) > 0)

datp <- normalize.percentControls(
  datp, 
  groupBy             = 'BARCODE',
  WELLTYPE_CODES 		  = c('LC'),     				      # WELLTYPE_CODE(S) of the control wells
  featureColumns      = cellcountfeatures_raw,		# columns to normalize, NULL: all numeric columns not in 'colnames.ignore'
  agg.function        = median,
  suffix              = '.prctCntr',
  colnames.ignore 	  = NULL
)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#  normalize all features: z-score relative to DMSO controls
#
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
datp <- normalize.zscoreControls_plateMean_globalPooledStdev(
  datp,
  groupBy        			    = c('expID'),					# these groups will be handled completely independently
  groupBy.mean        	  = c('BARCODE'),				# grouping for means (or medians)
  groupBy.pooled_sd 		  = NULL,   						# grouping for s.d. (or mad), pool the s.d. (or mad) into one value, NULL: around all plates
  WELLTYPE_CODES 			    = c('LC'),     	 			# WELLTYPE_CODE(S) of the control wells
  featureColumns      	  = features, 					# columns to normalize, NULL: all numeric columns not in 'colnames.ignore'
  suffix              	  = '.zGSLC',
  centerBeforeComputingSD = TRUE,							  # TRUE: subtract within each "groupBy.mean"-group the mean (or median), then compute s.d. // FALSE: compute s.d. on non-centered data
  robust              	  = FALSE,        			# FALSE: mean / sd, TRUE: median / mad
  colnames.ignore 		    = NULL
)
features_norm <- attr(datp, 'features.normalized')

# drop raw features from data frame, except cell count feature
features_raw_drop <- features[!grepl('CellCount_AllRetained', features)]
datp <- datp[, setdiff(names(datp), features_raw_drop)]

browse_table(data.frame.OverviewColumnStats.2(datp[, features_norm]), 'stats after normalization')
stopifnot(!any(is.na(datp$Treatment)))

#:::::::::::::::::::::::::::::::::::::::::::::::::
#
#  feature selection
#
#:::::::::::::::::::::::::::::::::::::::::::::::::
source('PhenoPrint/MasterThesis2017/featureSelection_mRMR.R')
columnReplicateID 	<- "Treatment"  # the column that contains the replicate information; all rows with the same value are considered replicates
cellcount_ft 		    <- 'CellCount_AllRetained'
min_cellcount       <- 100
nFeatures 			    <- 75 	# number of features to select (initially, the optimal number which is usually less than nFeatures will be determined below)
N_REPLICATES 		    <- 0    # > 0: randomly sample n replicates per treatment (to reduce influence of compounds with large number of replicates)

AUCeval.m 						          <- 10
AUCeval.nPairsReplicates 		    <- 0.2 
AUCeval.nPairsNonReplicates 	  <- 0.2
AUCeval.nPairsReplicatesMax 	  <- 10000
AUCeval.nPairsNonReplicatesMax 	<- 10000
AUCeval.distMethod 				      <- "pearson"
AUCeval.nThresholdsAUC 			    <- 100

DIR_OUT 						            <- 'C:/TEMP'

# -------------------------------
#  specify the control value
# -------------------------------
valueControl <- unique(datp[datp$WELLTYPE_CODE == 'LC', columnReplicateID])
stopifnot(length(valueControl) == 1)

# ----------
#  run mRMR
# ----------
FSres <- run_mRMR(
  dataset             = datp,
  features       	    = features_norm,    # all features to select from
  nFeatures 			    = nFeatures, 		    # number of features to be selected
  min_cellcount 	    = min_cellcount, 		# select features based on wells with at least this many cells
  cellcount_ft 		    = cellcount_ft, 
  N_REPLICATES 		    = N_REPLICATES,   	# > 0: randomly sample n replicates per treatment (to reduce influence of compounds with large number of replicates)
  columnReplicateID 	= columnReplicateID,
  valueControl        = valueControl,
  AUCeval.m 			                = AUCeval.m,
  AUCeval.nPairsReplicates 		    = AUCeval.nPairsReplicates, 
  AUCeval.nPairsNonReplicates 	  = AUCeval.nPairsNonReplicates,
  AUCeval.nPairsReplicatesMax 	  = AUCeval.nPairsReplicatesMax,
  AUCeval.nPairsNonReplicatesMax 	= AUCeval.nPairsNonReplicatesMax,
  AUCeval.distMethod 				      = AUCeval.distMethod,
  AUCeval.nThresholdsAUC 			    = AUCeval.nThresholdsAUC,
  DIR_OUT = DIR_OUT
)
features_norm_sel <- FSres$selectedFeatures[1:FSres$optN] # select the optimal number of features

#:::::::::::::::::::::::::::::::::::::::::::::::::
#
#  active calling based on selected feature set
#
#:::::::::::::::::::::::::::::::::::::::::::::::::
source('PhenoPrint/MasterThesis2017/active_calling.R')

stopifnot(!any(is.na(datp$Treatment)))

datp <- active_calling(datp, features_norm_sel)
hist(datp$Active_FractionOfRepl) # treatments with >= 50% active replicates are considered active

stopifnot(!any(is.na(datp$Treatment)))

#:::::::::::::::::::::::::::::::::::::::::::::::::
#
#  save preprocessed data to file
#
#:::::::::::::::::::::::::::::::::::::::::::::::::
attr(datp, 'features')           <- features_norm 
attr(datp, 'features_mRMR_all')  <- FSres$selectedFeatures
attr(datp, 'features_mRMR_optN') <- features_norm_sel

save(datp, file = fn_out)
cat('written to ',fn_out,'\n')

