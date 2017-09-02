# Project: Steffens-R-Projects
# 
# Author: SJAENSCH
###############################################################################

require(ggplot2)

source("E:/uhasselt/22yrsem/MThesis/scripts/R_code/stools.R")	

replicateCorrelation <- function(
		data.AllPlates, 		# well data
		features, 				# vector of feature names for which to calculate correlation
		groupBy = NULL,			# treat rows having different 'groupBy' value separately from each other
		replicateBy,  			# treat rows having the same 'replicateBy' value as replicates, thus replicates are defined grouping c(groupBy, replicateBy)
		axisBy 	= NULL,			# within replicates, these variables define which replicate goes on the x-axis and which on the y-axis; each axisBy value will be compared once against all other axisBy values
		shuffle = FALSE,        # if axisBy == NULL, shuffle the replicates to randomize which replicate goes on the x-axis and which on the y-axis
		summarizePerGroup 	= FALSE,
		summarizePerFeature = TRUE,
		fn.pdf              = NULL   # NULL: don't plot. NA or filename: x vs. y plots will be plotted to PDF file; if NA currently open pdf file will be used
)
{
	if(!is.null(fn.pdf) && !is.na(fn.pdf))
		fn.pdf <- init.newPDF(fn.pdf)
	
	if(length(axisBy) > 0)
		list.cor <- replicateCorrelation.byAxis(data.AllPlates, features, groupBy, replicateBy, axisBy, plot = !is.null(fn.pdf)) 
	else
		list.cor <- replicateCorrelation.allAgainstAll(data.AllPlates, features, groupBy, replicateBy, shuffle = shuffle, plot = !is.null(fn.pdf)) 
	
	if(!is.null(fn.pdf) && !is.na(fn.pdf))
		finish.pdf(fn.pdf, openPDF = TRUE)
	
	df.cor 	<- rbind.fill(list.cor)
	df.cor 	<- na.omit(df.cor)
	df.cor 	<- arrange(df.cor, -repl.cor)
	
	res 	<- list(df.cor = df.cor)
	
	# summarize each feature per group; pair count weighted average correlation
	if(summarizePerGroup)
	{	df.cor.agg <- ddply(
				df.cor, 
				c('feature', groupBy), 
				function(df) data.frame(
							repl.cor 		= weighted.mean(df$repl.cor, df$pair.cnt), 
							pair.cnt.avg 	= mean(df$pair.cnt), 
							n.groups 		= nrow(df))
		)
		df.cor.agg 	<- arrange(df.cor.agg, 	-repl.cor)
		res 		<- c(res, list(df.cor.groupSummary = df.cor.agg))
	}
	
	
	if(summarizePerFeature)
	{
		# summarize each feature; pair count weighted average correlation
		df.cor.agg2 <- ddply(
				df.cor, 
				'feature', 
				function(df) data.frame(
							repl.cor 		= weighted.mean(df$repl.cor, df$pair.cnt), 
							pair.cnt.avg 	= mean(df$pair.cnt), 
							n.groups 		= nrow(df))
		)
		df.cor.agg2 <- arrange(df.cor.agg2, 	-repl.cor)
		res 		<- c(res, list(df.cor.featureSummary = df.cor.agg2))
	}
	return(res)
}


replicateCorrelation.byAxis <- function(
		data.AllPlates, 		# well data
		features, 				# vector of feature names for which to calculate correlation
		groupBy = NULL,			# treat rows having different 'groupBy' value separately from each other
		replicateBy,  			# treat rows having the same 'replicateBy' value as replicates, thus replicates are defined grouping c(groupBy, replicateBy)
		axisBy,					# within replicates, these variables define which replicate goes on the x-axis and which on the y-axis; each axisBy value will be compared once against all other axisBy values
		plot  = FALSE   		# TRUE: x vs. y plots will be plotted to PDF file

)
{
	data.frame.assertColumnsExist(data.AllPlates, c(replicateBy, axisBy, groupBy, features))
	
	list.grps 			<- dlply(data.AllPlates[, c(replicateBy, axisBy, groupBy, features)], groupBy )
	col.axisByMerged 	<- c(paste(axisBy, '.x', sep = ''), paste(axisBy, '.y', sep = '')) # column names in the replicates-matched data.frame (see merge() call below)
	list.cor <- list()
	cnt <- 0
	for(i in seq_len(length(list.grps)))
	{
		df.i <- list.grps[[i]]
		list.df <- dlply(df.i, axisBy)
		
		
		if(length(list.df) <  2) {
			cat('No axis pairs in group:\n')
			str(df.i[1, groupBy, drop = FALSE])
			warning('No axis pairs in group (see above). Cannot compute correlation for this group.')
			next
		} 	
		
		for(a in 1 : (length(list.df) - 1)) # loop over 'x-axis' 
		{
			cat(a, ' / ', length(list.df), '\n', sep = '')
			da <- list.df[[a]]
			
			for(b in (a + 1) : length(list.df)) # loop over 'y-axis' 
			{
				db <- list.df[[b]]
				
				d <- merge(da, db, by = c(groupBy, replicateBy), suffixes = c('.x', '.y'))
				if(nrow(d) == 0)
					next
				
				repl.cor <- numeric(length(features))
				pair.cnt <- numeric(length(features))
				for(k in 1 : length(features)) # loop over feature
				{
					feature_k <- features[k]
					
					x <- d[, paste(feature_k, '.x', sep = '')]
					y <- d[, paste(feature_k, '.y', sep = '')]
					
					# explicitly remove NA values for accurate pair count
					{
						ix 	<- !is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y)
						x 	<- x[ix]
						y 	<- y[ix]
					}
					
					if(length(x) < 2)
						repl.cor[k] <- NA
					else					
						repl.cor[k] <- suppressWarnings(cor(x, y, use = 'all.obs')) # standard deviation might be zero ==> will result in NA correlation
					
					if(plot)
					{
						
						plotReplicates(x, y, feature_k, repl.cor[k], tlabel = '',
								xlabel = paste(t(da[1,axisBy]), collapse = ' // '),
								ylabel = paste(t(db[1,axisBy]), collapse = ' // ')
						)
					}
					
					pair.cnt[k] <- length(x)
				}
				
				suppressWarnings( # short variable row names will be discarded
						d.cor <- cbind(
								data.frame(feature = features, repl.cor = repl.cor, pair.cnt = pair.cnt, stringsAsFactors = FALSE),
								d[1, c(groupBy, col.axisByMerged)]
						)
				)
				cnt <- cnt + 1
				list.cor[[cnt]] <- d.cor
			}
		}
	}
	
	return(list.cor)
}	


replicateCorrelation.allAgainstAll <- function(
		data.AllPlates, 		# well data
		features, 				# vector of feature names for which to calculate correlation
		groupBy = NULL,			# treat rows having different 'groupBy' value separately from each other
		replicateBy,  			# treat rows having the same 'replicateBy' value as replicates, thus replicates are defined grouping c(groupBy, replicateBy)
		shuffle = FALSE,		# shuffle the replicates to randomize which replicate goes on the x-axis and which on the y-axis
		plot  	= FALSE   		# if TRUE: x vs. y plots will be plotted to PDF file
)
{
	list.grps	<- dlply(data.AllPlates[, c(replicateBy, groupBy, features)], groupBy )
	list.cor 	<- list()
	cnt 		<- 0
	for(i in seq_len(length(list.grps))) # loop over groups; they will be treated completely separately from each other
	{
		cat('Group ', i, ' / ', length(list.grps), '\n', sep = '')
		df.i <- list.grps[[i]]
		
		df.i$ROWIDX__ <- 1 : nrow(df.i) # assign row index
		
		list.repl <- dlply(df.i, replicateBy) # group by replicates; each element in the list is a data.frame of replicates
		
		#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
		#
		#     get the row indices of the replicates to be compared
		#     - compare each replicate once against all other replicates
		#
		#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
		{
			x.rowidx <- integer(10000)
			y.rowidx <- integer(10000)
			cnt2 	 <- 0
			for(j in 1 : length(list.repl)) # loop over replicates 
			{
				df.j <- list.repl[[j]]
				if(nrow(df.j) < 2)
					next
				
				if(shuffle)
					df.j <- df.j[sample(1:nrow(df.j), size = nrow(df.j), replace = FALSE),]
				
				for(a in 1 : (nrow(df.j) - 1)) # loop over 'x-axis' 
				{
					for(b in (a + 1) : nrow(df.j)) # loop over 'y-axis' 
					{
						cnt2 <- cnt2 + 1
						x.rowidx[cnt2] <- df.j$ROWIDX__[a] 
						y.rowidx[cnt2] <- df.j$ROWIDX__[b] 
					}
				}
			}
			x.rowidx <- x.rowidx[1:cnt2]
			y.rowidx <- y.rowidx[1:cnt2]
		}
		
		# compute correlation for each feature
		{
			repl.cor <- numeric(length(features))
			pair.cnt <- numeric(length(features))
			for(k in 1 : length(features)) # loop over feature
			{
				feature_k <- features[k]
				
				x <- df.i[x.rowidx, feature_k]
				y <- df.i[y.rowidx, feature_k]
				
				# explicitly remove NA values for accurate pair count
				{
					ix 	<- !is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y)
					x 	<- x[ix]
					y 	<- y[ix]
				}
				
				pair.cnt[k] <- length(x)
				if(length(x) < 2)
					repl.cor[k] <- NA
				else
				{
					repl.cor[k] <- cor(x, y, use = 'all.obs' )
					
					if(plot)
					{
						if(length(groupBy) > 0)
							tlabel <- paste(t(df.i[1, groupBy]), collapse = ' // ')
						else
							tlabel <- ''
						
						plotReplicates(x, y, feature_k, repl.cor[k], tlabel = tlabel)
						
					}
				}
			}
		}	
		suppressWarnings( # short variable row names will be discarded
				d.cor <- cbind(
						data.frame(feature = features, repl.cor = repl.cor, pair.cnt = pair.cnt, stringsAsFactors = FALSE),
						df.i[1, c(groupBy), drop = FALSE]
				)
		)
		cnt <- cnt + 1
		list.cor[[cnt]] <- d.cor
	}
	
	return(list.cor)
}

plotReplicates <- function(x, y, featureName, correlation, tlabel = '', xlabel = featureName, ylabel = featureName)
{
	lim 	<- range(c(x,y))
	margin 	<- max(abs(lim)) * 0.05
	lim 	<- lim + c(-margin, +margin)
	
	dgg <- data.frame(x = x, y = y)
	sp <- ggplot(dgg, aes(x = x, y = y)) + 
			geom_hex(bins = 75) + 
			scale_fill_gradientn(colours = rainbow(7)) +
			geom_abline(intercept = 0, slope = 1) + 
			ggtitle(paste(featureName, '\nCorrelation = ', sprintf('%.4f', correlation), '\n', tlabel)) +
			xlab(xlabel) +
			ylab(ylabel) +
			xlim(lim) +
			ylim(lim) +
			#coord_cartesian(xlim = lim, ylim = lim) +
			theme(aspect.ratio = 1) # square plotting area
	print(sp)
}
##::::::::::::::::::::::::::::::::::::::
##
##   compute replicate variance
##
##::::::::::::::::::::::::::::::::::::::
#repl.var 	<- numeric(length(features))
#for(k in 1 : length(features))
#{
#	feature_k <- features[k]
#	cat(k, ' / ', length(features), ': ', feature_k, '\n', sep = '')
#	
#	repl.var_i 	<- numeric(100000)
#	repl.cnt_i 	<- numeric(100000)
#	cnt <- 0
#	for(i in 1 : length(list.repl))
#	{
#		df.i <- list.repl[[i]]
#		
#		values <- df.i[, feature_k] 
#		values <- values[!is.na(values)]
#
#		if(length(values) < 2)
#			next
#		
#		cnt <- cnt + 1
#		repl.var_i[cnt] <- var(values)
#		repl.cnt_i[cnt] <- length(values)
#		
#	}
#	repl.var_i <- repl.var_i[1:cnt]
#	repl.cnt_i <- repl.cnt_i[1:cnt]
#	
#	repl.var[k] <- sum(repl.cnt_i * repl.var_i) / sum(repl.cnt_i)  # replicate count weighted average variance
#}
