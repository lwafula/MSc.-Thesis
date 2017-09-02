#LWAFULA (07072017)
#the feature plots scatter plots have been commented out to speed up, 


library(MRMRFS)

source("E:/uhasselt/22yrsem/MThesis/scripts/R_code/stools.R")	
source("E:/uhasselt/22yrsem/MThesis/scripts/R_code/ReplicateCorrelationVariance.R")	

removeRows <- function(df, ix, msg)
{
	if(is.logical(ix))
		ix <- which(ix)
	
	if(length(ix) > 0)
	{
		n0 <- nrow(df)
		df <- df[-ix, ]
		cat('Removed ', length(ix), ' / ',n0,' rows: ', msg, '\n', sep = '')
	}
	return(df)
}


run_mRMR <- function(
		dataset,
		features,				        # all features to select from
		nFeatures, 				      # number of features to be selected
		min_cellcount = 100, 	  #  select features based on wells with at least this many cells
		cellcount_ft, 
		N_REPLICATES,   		    # > 0: randomly sample n replicates per treatment (to reduce influence of compounds with large number of replicates)
		columnReplicateID = "Treatment",
		valueControl,			      # value in column 'columnReplicateID' for the controls, must be unique for all controls
		AUCeval.m,
		AUCeval.nPairsReplicates, 
		AUCeval.nPairsNonReplicates,
		AUCeval.nPairsReplicatesMax,
		AUCeval.nPairsNonReplicatesMax,
		AUCeval.distMethod 				= 'Pearson',
		AUCeval.nThresholdsAUC 		= 100,
		
		label 		= '',
		DIR_OUT 	= 'C:/TEMP' # MRMR_<label> will be appended as sub-folder

)
{
	
	stopifnot(length(valueControl) == 1)
	
	#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
	#
	#   remove wells with low cell count
	#
	#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
	cc <- dataset[, cellcount_ft]
	datp <- removeRows(dataset, cc < min_cellcount, msg = paste0('less than ',min_cellcount,' cells'))
	
	#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	#
	# 	remove rows with NA (selectActiveIds does not accept NAs)
	#
	#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	ix.complete <- complete.cases(dataset[, c(columnReplicateID, features)]) # select non-NA rows only based on the columns required for MRMR
	dataset <- removeRows(dataset, !ix.complete, msg = paste0("containing NA's"))
	
	#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	#
	# 	select required columns
	#
	#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	cols <- NULL
	dataset <- dataset[, c(columnReplicateID, cols, features, 'WELLTYPE_CODE')]  
	
	#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	#
	#      select max N replicates per cmpd/concentration
	#      (random sampling)
	#
	#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	if(N_REPLICATES > 0)
	{
		cat('Selecting ',N_REPLICATES, ' random replicates per treatment.\n')
		n0 <- nrow(dataset)
		dataset <- data.frame.samplePerGroup(
				dataset, 
				groupBy = columnReplicateID, 
				size 	= N_REPLICATES, 
				replace = FALSE, 
				errorOn.notEnoughSamplesInGroup = FALSE
		)
		n1 <- nrow(dataset)
		cat(n1, ' / ', n0, ' rows left\n')	
	}
	
	#:::::::::::::::::::::::::::::::::::::::::::::::::::
	#
	#  remove treatments with only 1 replicate
	#
	#:::::::::::::::::::::::::::::::::::::::::::::::::::
	replCount 	<- table(dataset[, columnReplicateID])
	singletons 	<- replCount[replCount == 1]
	dataset 	<- removeRows(dataset, dataset[, columnReplicateID] %in% names(singletons), msg = paste0("treatments with only 1 replicate."))
	
	cat('Number of treatments with more than 1 replicate: ',length(unique(dataset[, columnReplicateID])),'\n')
	
	DIR_OUT <- file.path(DIR_OUT, paste0('mRMR',label))
	ensure.dirExists(DIR_OUT)
	fn.pdf <- init.newPDF(file.makeFilenameUnique(file.path(DIR_OUT, paste0('MRMR',label,'.pdf')), appendIfAlreadyUnique = TRUE))
	
	#:::::::::::::::::::::::::::::::::::::::::::::::::::
	#
	#    select active treatments
	#
	#:::::::::::::::::::::::::::::::::::::::::::::::::::
	tic()
	
	ptm <- proc.time()
	n0 <- nrow(dataset)
	activeSelection <- selectActiveIds(dataset, cols, features, columnReplicateID, columnReplicateID,  valueControl)
	
	proc.time() - ptm
	
	# the data with active rows 
	dataset_active <- activeSelection$hcsData # note: "dataset_active" is a data.table, while "dataset" is a data.frame
	n1 <- nrow(dataset_active)
	
	cat(n1, ' / ', n0, ' rows remaining after active selection.\n')
	cat('Number of active treatments with more than 1 replicate: ', activeSelection$nIds,'\n')
	
	
	# tmp <- data.frame(dataset_active)
	# debug.write.csv(counttable(tmp[, columnReplicateID]), name = paste0('active_treatments ',label))
	
	#:::::::::::::::::::::::::::::::::::::::::::::::::::
	#
	#    MRMR computation
	#
	#:::::::::::::::::::::::::::::::::::::::::::::::::::
	message("MRMR method computation ...")
	ptm <- proc.time()
	MRMR_output <- MRMR(dataset_active, features, columnReplicateID, nFeatures)
	proc.time() - ptm
	
	selectedFeatures <- MRMR_output$selectedFeatures
	
	# relevance (V)
	{
		fnOut_V <- file.path(DIR_OUT, paste0('MRMR',label,'_V.csv'))
		write.csv(MRMR_output$stat, file = fnOut_V, row.names = FALSE)
		
		V <- sort(MRMR_output$stat$V, decreasing = TRUE)
		plot(V, xlab = 'Feature Index', ylab = 'Relevanz (V)', main = 'Relevanz',pch=19,cex=0.75)
	}
	
	# selectionStatistic
	{
		fnOut_selStat <- file.path(DIR_OUT, paste0('MRMR',label,'_selectionStatistic.csv'))
		d0 <- data.frame(
				rank = 1 : length(selectedFeatures), 
				feature = selectedFeatures, 
				selectionStatistic = MRMR_output$selectionStatistic
		)
		write.csv(d0, file = fnOut_selStat, row.names = FALSE)
		
		plot(MRMR_output$selectionStatistic, 
				xlab = 'Selected Feature Index', 
				ylab = 'MRMR Criterion (V = relevance - (relevance * redundancy))', 
				main = 'MRMR Criterion (V = relevance - (relevance * redundancy))',pch=19,cex=0.65)
	}
	
	#:::::::::::::::::::::::::::::::::::::::::::::::::::
	#
	#    correlation between features
	#
	#:::::::::::::::::::::::::::::::::::::::::::::::::::
	
	computationCor <- computeCor(dataset_active, selectedFeatures = selectedFeatures, nFeatures)
	
	plotCor(computationCor$maxCorAll, type = c("maxCorAll"))

	plotCor(computationCor$maxCorPrevious, type = c("maxCorPrevious"))

	df_featureCor <- data.frame(
			rank	= 1 : length(selectedFeatures),
			feature = selectedFeatures, 
			maxCorToPrevFeature = c(0, computationCor$maxCorPrevious), 
			maxCorToAllFeatures = computationCor$maxCorAll
	)
	
	#:::::::::::::::::::::::::::::::::::::::::::::::::::
	#
	#    How many features are sufficient?
	#
	#:::::::::::::::::::::::::::::::::::::::::::::::::::
	
	# ROC-AUC approach
	{
		tic()
		#nFeaturesConsidered <- unique(round(10^(seq(log10(2), log10(length(features)), length.out = 20))))
		nFeaturesConsidered <- 2 : length(selectedFeatures)
		
		optimalNFeatAUCPearson <- evaluationAUC(
				datasetActive 			    = as.data.frame(dataset_active),
				features 				        = selectedFeatures,
				nFeaturesConsidered 	  = nFeaturesConsidered,
				columnReplicateID 		  = columnReplicateID,
				m 						          = AUCeval.m,
				nPairsReplicates 		    = AUCeval.nPairsReplicates, 
				nPairsNonReplicates 	  = AUCeval.nPairsNonReplicates,
				nPairsReplicatesMax 	  = AUCeval.nPairsReplicatesMax,
				nPairsNonReplicatesMax 	= AUCeval.nPairsNonReplicatesMax,
				distMethod 				      = AUCeval.distMethod,
				nThresholdsAUC 			    = AUCeval.nThresholdsAUC
		)
		plot(x = optimalNFeatAUCPearson)
		optN <- optimalNFeatAUCPearson$optimalNFeatures # optimal number of features
		toc()
		cat('AUC-based evaluation of number of features finished\n')
	}
	
	#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	#
	# print formatted feature vector (as R-code) to console
	#
	#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	cat( '# optimal number of features: ', optN, '\n# ',label,'\n\n',
			paste(
					"'",
					selectedFeatures, 
					c(rep("', #", times = length(selectedFeatures)-1), "' #"),
					1:length(selectedFeatures),
					sep = '', 
					collapse = '\n'
			), 
			'\n', 
			sep = ''
	)
	
	#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	#
	#   scatterplots of replicates for the selected features
	#   (for visualization purposes only; maybe commented out for speed-up)
	#
	#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	save(dataset_active, selectedFeatures, columnReplicateID, file = 'C:/TEMP/mrmr.Rdata')
   replCor <- replicateCorrelation(
   		as.data.frame(dataset_active), 		# well data
   		selectedFeatures, 					      # vector of feature names for which to calculate correlation
   		groupBy 	= NULL, 				        # treat rows having different 'groupBy' value separately from each other
   		replicateBy = columnReplicateID,	# treat rows having the same 'replicateBy' value as replicates, thus replicates are defined grouping c(groupBy, replicateBy)
   		axisBy 		= NULL,					        # within replicates, these variables define which replicate goes on the x-axis and which on the y-axis; each axisBy value will be compared once against all other axisBy values
   		                                  # NULL: no axis; match each replicate against each other replicate
   		summarizePerGroup 	= TRUE,
   		summarizePerFeature = TRUE,
   		fn.pdf              = NA			# NA: use currently open pdf
   )
 	
 	finish.pdf(fn.pdf, openPDF = TRUE)
 	
	
	#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	#
	#   return results
	#
	#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	r <- list(
			label 			     = label, 
			selectedFeatures = selectedFeatures, 
			allFeatures      = features,
			optN 			       = optN
	)
	
	return(r)
}

