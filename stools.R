require(matlab, quietly = TRUE)
require(plyr)
require(data.table)

#:::::::::::::::::::::::::::::::::::::::::::
#
#  file tools
#
#:::::::::::::::::::::::::::::::::::::::::::

ensure.dirExists <- function(path)
{
  fparts <- fileparts(path)
  isdir  <- file.info(path)$isdir
  
  if( (is.na(isdir) && fparts$ext == "" ) || (!is.na(isdir) && isdir) )
    path_ <- path 
  else
    path_ <- fparts$pathstr
  
  # remove trailing slash, file.exists will fail otherwise on windows!!!
  path_ <- gsub('/$','',path_)
  #cat('path_: ',path_,'\n')
  if(!file.exists(path_))
    dir.create(path_,recursive=TRUE)
  
  return(path_)
}

ensure.emptyDirExists <- function(path)
{
  fparts <- fileparts(path)
  isdir  <- file.info(path)$isdir
  
  if(isdir || (is.na(isdir) && fparts$ext == "" ))
    path_ <- path 
  else
    path_ <- fparts$pathstr
  
  
  if(file.exists(path_))
    unlink(path_,recursive = TRUE)
  
  dir.create(path_,recursive=TRUE)
  return(path)
}


file.makeFilenameValid <- function(name)
{
  return(gsub('\\/|\\\\|:|\\*|\\?|"|<|>|\\|', '_', name))
}

file.makeFilenameUnique <- function(fn, sep = '_', startAt = 1, appendIfAlreadyUnique = FALSE)
{
  fp 	<- fileparts(fn)
  idx <- startAt
  
  if(appendIfAlreadyUnique)
    fn 	<- file.path(fp$pathstr, paste(fp$name, sep, idx, fp$ext, sep = ''))
  
  while(file.exists(fn))
  {
    fn 	<- file.path(fp$pathstr, paste(fp$name, sep, idx, fp$ext, sep = ''))
    idx <- idx + 1
  }
  return(fn)
}

# removes the file if exists
# if removing fails, "file.makeFilenameUnique" will be called and a new file name returned
file.clear <- function(fn, ...)
{
  if(file.exists(fn))
  {
    succeeded <- FALSE
    tryCatch(
      {
        suppressWarnings( file.remove(fn))				
        succeeded <- !file.exists(fn)
      },
      error = function(e){}
    )
    if(!succeeded)
    {
      fn <- file.makeFilenameUnique(fn)
    }
  }
  stopifnot(!file.exists(fn))
  return(fn)
}

#:::::::::::::::::::::::::::::::::::::::::::
#
#  debugging tools
#
#:::::::::::::::::::::::::::::::::::::::::::

debug.write.csv <- function(
  df, 
  file 		= NULL, 
  name 		= NULL, # for convenience: construct file using this name and the path from the specified (or default) file 
  openFile 	= TRUE, 
  row.names 	= FALSE,
  na 			= "NA",
  quote 		= TRUE,
  default_file = 'C:/TEMP/aTable.csv'
)
{
  if(is.null(file))
  {
    file <- file.makeFilenameUnique(default_file)	
    if(is.null(name))
      name <- deparse(substitute(df)) # for convenience: construct file using this name and the path from the specified (or default) file 
  } else {
    name <- NULL
  }
  
  if(!is.null(name))
  {
    succeeded <- FALSE
    tryCatch(
      {
        name <- gsub('\\/|\\\\|:|\\*|\\?|"|<|>|\\|', '_', name) # remove characters not allowed in windows filename (replace by _)		
        file <- file.makeFilenameUnique(file.path(file.getPath(file), paste(name[1], '.csv', sep = '')))
        
        succeeded <- TRUE
      },
      error = function(e){}
    )
    if(!succeeded)
    {
      file <- file.makeFilenameUnique(default_file)
    }
    
  }	
  succeeded <- FALSE
  tryCatch(
    {
      write.csv(df, file = file, row.names = row.names, na = na, quote = quote  )
      succeeded <- TRUE
    },
    error = function(e){}
  )
  
  if(!succeeded)
  {
    tryCatch(
      {
        file <- file.makeFilenameUnique(file)
        write.csv(df, file = file, row.names = row.names, na = na, quote = quote )
        succeeded <- TRUE
      },
      error = function(e){}
    )
  }
  
  if(!succeeded)
  {
    file <- file.makeFilenameUnique(default_file)
    write.csv(df, file = file, row.names = row.names, na = na, quote = quote )
  }
  
  
  if(openFile)
    browseURL(file)
  return(file)
}

browse_table <- function(df, title = 'aTable')
{
  # debug.write.csv(df, name = title)
  View(df, title = title)
}

#:::::::::::::::::::::::::::::::::::::::::::
#
#  data analysis & manupulation tools
#
#:::::::::::::::::::::::::::::::::::::::::::
data.frame.OverviewColumnStats.2 <- function(df)
{
  df <- as.data.frame(df)
  
  cn <- colnames(df)
  list.df.stats <- list()
  for (i in 1:length(cn))
  {
    cn_i <- cn[i]
    x <- df[, cn_i]
    n <- length(x)
    if (is.numeric(x))
    {
      isNum 	<- TRUE
      nix.na 	<- which(is.na(x))
      nix.NaN <- which(is.nan(x))
      nix.Inf <- which(is.infinite(x))
      nix.0 	<- which(x == 0)
      x 		<- x[!is.na(x) & is.finite(x)]
      v 	  <- var(x)
      cv 		<- sd(x) / mean(x)
      mad_	<- mad(x)
      mean_	<- mean(x)
      median_	<- median(x)
    } else{
      isNum 	<- FALSE
      nix.na 	<- which(is.na(x))
      nix.NaN <- integer(0)
      nix.Inf <- integer(0)
      nix.0 	<- integer(0)
      v 	  	<- NA
      cv 		<- NA
      mad_	<- NA
      mean_	<- NA
      median_	<- NA
    }
    
    list.df.stats[[i]] <- data.frame(
      colnames = cn_i,
      isNumeric 	= isNum,
      n 			= n,
      n_NA_Inf_0  = length(nix.na) + length(nix.Inf) + length(nix.0),
      n_NA 		= length(nix.na),
      n_NaN 		= length(nix.NaN),
      n_Inf 		= length(nix.Inf),
      n_0 		= length(nix.0),
      variance   	= v,
      CoV			= cv,
      mad 		= mad_,
      mean_		= mean_,
      median_		= median_
    )
  }
  df.stats <- rbind.fill(list.df.stats)
  df.stats <- arrange(df.stats, n_NA_Inf_0)
  return(df.stats)
}

data.frame.assertColumnsExist <- function(df, col.names, error = TRUE)
{
  if(length(col.names) == 0)
    return(NULL)
  
  cn 			<- colnames(df)
  col.names 	<- unique(col.names)
  missing 	<- col.names[!(col.names %in% cn)]
  
  if(length(missing) > 0)
  {
    name.df <-  deparse(substitute(df))
    
    cat(sprintf('The following columns are not in the data.frame "%s":\n%s\n', name.df, paste('\t','"', missing, '"', sep = '', collapse = '\n')))
    if(error)
      stop(sprintf('The following columns are not in the data.frame "%s" (see full list printed above):\n%s\n', name.df, paste('\t','"', missing, '"', sep = '', collapse = '\n')))
    else
      warning(sprintf('The following columns are not in the data.frame "%s" (see full list printed above):\n%s\n', name.df, paste('\t','"', missing, '"', sep = '', collapse = '\n')))
  }
  return(missing)
}


aggregate.dataframe.2 <- function(
  df,                                 			 	# the data 
  unique.var  	= character(0),						 	# store the first value of these variables within each group along with the aggregated values; this is equivalent to applying an aggregation function that returns the first item of the input vector
  agg.var     	= NULL,    							 	# aggregate these variables (i.e., columns), NULL: all numeric columns not in c(unique.var, groupBy)
  groupBy     	= unique.var,           			 	# group data.frame by these variables ==> aggregate each group
  agg.func    	= c('mean',          'sd'),        	 	# apply these aggregation functions to each column specified in agg.var        
  agg.colname 	= c('mean_<agg.var>','sd_<agg.var>') 	# name of the output columns, with <agg.var> being the place holder for the name of the input column
)
{
  if(is.null(df))
    return(NULL)
  
  unique.var <- unique(c(groupBy, unique.var))
  
  if(is.null(agg.var))
  {
    cn 		<- colnames(df)
    agg.var <- cn[as.logical(lapply (df, is.numeric )) & !(cn %in% c(unique.var, groupBy))]
  }
  
  data.frame.assertColumnsExist(df, agg.var)
  
  require(plyr)
  agg <- function(df, unique.var, agg.var, agg.func, agg.colname)
  {
    if(length(agg.func) != length(agg.colname))
      stop('agg.func and agg.colname must have same length')
    
    #..................................................................
    #
    #  get values for columns with unique values
    #  (use value in first row instead of actually calling 'unique')
    #
    #..................................................................
    df.agg <- df[1,unique.var, drop = FALSE]
    
    #..................................................................
    #
    #   apply the aggregation function(s)
    #
    #..................................................................
    for(j in 1 : length(agg.func))
    {
      #..............................................
      #
      #  instantiate names for agg columns
      #
      #..............................................
      col.names <- character(length(agg.var))
      for(k in 1 : length(agg.var))
        col.names[k] <- gsub('<agg.var>',agg.var[k], agg.colname[j])
      
      #..................................................................
      #
      #  compute aggregation values
      #
      #..................................................................
      
      n <- length(agg.var)
      a <- numeric(n)
      for(i in 1 : n)
      {
        x 		<- df[,agg.var[i]]
        x 		<- x[!is.infinite(x) & !is.na(x)]
        a[i] 	<- do.call(agg.func[j], list(x))
      }
      df.tmp			 <- as.data.frame(t(a), stringsAsFactors = FALSE)
      colnames(df.tmp) <- col.names 
      df.agg 			 <- cbind.data.frame( df.agg, df.tmp)
    }		
    return(df.agg)
  }
  
  list.df <- dlply (df, groupBy, .fun = agg, unique.var, agg.var, agg.func, agg.colname)
  df.agg  <- rbind.fill(list.df)
  if(!is.null(df.agg))
    attr(df.agg, 'agg.func') <- agg.func
  
  return(df.agg)
}


# allows to specify multiple groups of variables to be aggregated with each its own aggregation function(s)
aggregate.dataframe.3 <- function(
  df,                                 		# the data 
  unique.var  = character(0),					# store these variables (unique within each group!) along with the aggregation values
  groupBy     = unique.var,           		# group data.frame by these variables ==> aggregate each group
  list.agg    = list( 						# list of list(var = myvariables, func = myfunctions, colname = mycolnames), var == NULL ==> use all numeric columns 
    list(                        	
      var     = NULL,    								# aggregate these variables (i.e., columns), NULL: all numeric columns not in c(unique.var, groupBy)
      func    = c('mean',          'sd'),             # aggregation functions to be applied (given as strings, not as function handles)    
      colname = c('mean_<agg.var>','sd_<agg.var>')    # naming patterns of the output columns, <agg.var> = placeholder for the resp. input column
    )
  )
)
{
  if(is.null(df))
    return(NULL)
  
  unique.var <- unique(c(groupBy, unique.var))
  
  #::::::::::::::::::::::::::::::::::::::::::::::
  #
  #  1. check input
  #  2. agg.var NULL ==> find all numeric variables
  #
  #::::::::::::::::::::::::::::::::::::::::::::::
  for(k in 1 : length(list.agg)) # loop over groups of variables to be aggregated with the same function(s)
  {
    s <- list.agg[[k]]
    agg.var 	<- s$var
    agg.func	<- s$func
    agg.colname	<- s$colname
    
    if(length(agg.func) != length(agg.colname))
      stop(sprintf('agg.func and agg.colname must have same length (%d. entry in list.agg)', k))
    
    data.frame.assertColumnsExist(df, agg.var)
    
    if(is.null(agg.var))
    {
      cn 		<- colnames(df)
      s$var 	<- cn[as.logical(lapply (df, is.numeric )) & !(cn %in% c(unique.var, groupBy))]
    }
    list.agg[[k]] <- s
  }
  
  #::::::::::::::::::::::::::::::::::::::::::::::
  #
  #  aggregation of a single group
  #
  #::::::::::::::::::::::::::::::::::::::::::::::
  require(plyr)
  agg <- function(df, unique.var, list.agg)
  {
    #..................................................................
    #
    #  get values for columns with unique values
    #  (use value in first row instead of actually calling 'unique')
    #
    #..................................................................
    df.agg <- df[1,unique.var, drop = FALSE]
    
    #..................................................................
    #
    #  apply the aggregation function(s)
    #
    #..................................................................
    for(k in 1 : length(list.agg)) # loop over groups of variables to be aggregated with the same function(s)
    {
      s <- list.agg[[k]]
      agg.var 	<- s$var
      agg.func	<- s$func
      agg.colname	<- s$colname
      
      
      for(j in 1 : length(agg.func)) # loop over aggregation functions to be applied on the current variable group
      {
        #..............................................
        #
        #  instantiate names for agg columns
        #
        #..............................................
        col.names <- character(length(agg.var))
        for(p in 1 : length(agg.var))
          col.names[p] <- gsub('<agg.var>',agg.var[p], agg.colname[j])
        
        #..................................................................
        #
        #  compute aggregation values
        #
        #..................................................................
        
        n <- length(agg.var)
        a <- numeric(n)
        for(i in 1 : n)
        {
          x 		<- df[,agg.var[i]]
          x 		<- x[!is.infinite(x) & !is.na(x)]
          a[i] 	<- do.call(agg.func[j], list(x))
        }
        df.tmp			 <- as.data.frame(t(a), stringsAsFactors = FALSE)
        colnames(df.tmp) <- col.names 
        df.agg 			 <- cbind.data.frame( df.agg, df.tmp)
      }
    }
    return(df.agg)
  }
  
  list.df <- dlply(df, groupBy, .fun = agg, unique.var, list.agg)
  df.agg  <- rbind.fill(list.df)
  
  if(!is.null(df.agg))
    attr(df.agg, 'agg.func') <- agg.func
  
  return(df.agg)
}

data.frame.setAttributes.noOverwrite <- function(df, myAttributes)
{
  attr.df 		<- attributes(df)
  myNames 		<- names(myAttributes)
  dfNames 		<- names(attr.df)
  attributes(df) 	<- c(attr.df, myAttributes[!(myNames %in% dfNames)]) # df's attributes + myAttributes that are not in df
  return(df)
}

data.frame.assertColumnsExist <- function(df, col.names, error = TRUE)
{
  if(length(col.names) == 0)
    return(NULL)
  
  cn 			<- colnames(df)
  col.names 	<- unique(col.names)
  missing 	<- col.names[!(col.names %in% cn)]
  
  if(length(missing) > 0)
  {
    name.df <-  deparse(substitute(df))
    
    cat(sprintf('The following columns are not in the data.frame "%s":\n%s\n', name.df, paste('\t','"', missing, '"', sep = '', collapse = '\n')))
    if(error)
      stop(sprintf('The following columns are not in the data.frame "%s" (see full list printed above):\n%s\n', name.df, paste('\t','"', missing, '"', sep = '', collapse = '\n')))
    else
      warning(sprintf('The following columns are not in the data.frame "%s" (see full list printed above):\n%s\n', name.df, paste('\t','"', missing, '"', sep = '', collapse = '\n')))
  }
  return(missing)
}

data.frame.getColumns <- function(
  df, 
  col,			      # names of columns to be selected 
  error = TRUE   	# TRUE: throw error if there are missing columns, FALSE: drop missing columns (and warning)
)
{
  col <- unique(col)
  missing <- data.frame.assertColumnsExist(df, col, error = error)
  if(length(missing) > 0)
    col <- col[!(col %in% missing)]
  return(df[, col])
}

counttable <- function(x, colnameOut.count = NULL, rm.zeroCounts = TRUE)
{
  df <- data.frame(table(x))
  if(rm.zeroCounts)
    df <- df[df$Freq > 0,]
  
  # rename 'Freq' to user-specified name
  if(!is.null(colnameOut.count))
    colnames(df)[colnames(df) == 'Freq'] <- colnameOut.count
  
  return(df)
}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#   data normalization methods
#
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               normalize (generic)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
normalize.generic.singleGroup <- function(
  df.welldata, 
  colname.WELLTYPE_CODE 	= 'WELLTYPE_CODE',								# name of the column that takes the role of 'well type' (to distinugish controls from samples)
  WELLTYPE_CODES 			= c('LC'),      								# WELLTYPE_CODE(S) of the control wells
  featureColumns      	= NULL, 										# columns to normalize, NULL: all numeric columns not in 'colnames.ignore'
  func  					= function(x, cntr) {100 * x / median(cntr)},	# normalize using this function
  suffix              	= '.prctCntr',
  colnames.ignore 		= c("EXPERIMENT_ID","EXPERIMENT_NAME","PLATE_ID","BARCODE","SEQUENCE_IN_RUN","PLATE_INFO","VALIDATE_STATUS",
                        "APPROVE_STATUS","REMARKS","ROW_NR","COL_NR","WELLTYPE_CODE","COMPOUND_TY","COMPOUND_NR",
                        "CONCENTRATION","IS_VALID","FEATURE_NAME","FEATURE_ALIAS", colname.WELLTYPE_CODE),
  stopIfNoCntrFound = TRUE
)
{
  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  #
  #   get controls wells
  #
  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ix 			<- df.welldata[, colname.WELLTYPE_CODE] %in% WELLTYPE_CODES
  df.cntr 	<- df.welldata[ix, ]
  if(nrow(df.cntr) == 0)
  {
    if(stopIfNoCntrFound)
    {
      print(df.welldata)
      debug.write.csv(df.welldata, name = 'No_controls')
      stop(paste0('No controls of WELLTYPE_CODE ',paste0(WELLTYPE_CODES, collapse = '/'),' for normalization found for data printed above. see also excel table.'))
    } else {
      print(df.welldata[1,])
      warning(paste0('No controls of WELLTYPE_CODE ',paste0(WELLTYPE_CODES, collapse = '/'),' for normalization found for data printed above (1st row).'))
      return(df.welldata)
    }
  }
  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  #
  #  normalize the feature columns (samples and controls)
  #
  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  if(is.null(featureColumns))
    cnames <- colnames(df.welldata)	
  else
    cnames <- as.character(featureColumns)
  
  for(i in 1:length(cnames))
  {
    colname.i <- cnames[i]
    if(is.null(featureColumns) && colname.i %in% colnames.ignore)
      next
    
    x <- df.welldata[,colname.i]
    if(!is.numeric(x))
      next
    
    colname.i_norm 				 <- paste(colname.i, suffix, sep='')
    
    cntr 						 <- df.cntr[,colname.i]
    cntr 					 	 <- cntr[!is.na(cntr) & is.finite(cntr)]
    
    if(length(cntr) == 0)
      stop(sprintf('All controls value are NA or infinite in column %s', colname.i))
    
    x.norm 						 <- func(x, cntr)
    
    if(any(is.infinite(x.norm)))
    {
      cat(colname.i,'\n')
      cat('cntr:   ',  paste(cntr, 	collapse = ',\t'),'\n')
      cat('x:      ',  paste(x,     	collapse = ',\t'),'\n')
      cat('x.norm: ',  paste(x.norm, 	collapse = ',\t'),'\n')
      cat('------------------------------------------------\n')
      warning(sprintf('Infinite value produced when normalizing column %s against %s wells. Replacing by NaN\n', colname.i, paste0(WELLTYPE_CODES, collapse = '/')))
      x.norm[is.infinite(x.norm)] <- NaN
    }
    df.welldata[,colname.i_norm] <- x.norm
  }
  return(df.welldata)
}

normalize.generic <- function(
  df.welldata, 
  groupBy             	= c('BARCODE'),
  colname.WELLTYPE_CODE 	= 'WELLTYPE_CODE',								# name of the column that takes the role of 'well type' (to distinugish controls from samples)
  WELLTYPE_CODES 			= c('LC'),      # WELLTYPE_CODE(S) of the control wells
  featureColumns      	= NULL, 		# columns to normalize, NULL: all numeric columns not in 'colnames.ignore'
  func  					= function(x, cntr) {100 * x / median(cntr)},	# normalize using this function
  suffix              	= '.prctCntr',
  colnames.ignore 		= c("EXPERIMENT_ID","EXPERIMENT_NAME","PLATE_ID","BARCODE","SEQUENCE_IN_RUN","PLATE_INFO","VALIDATE_STATUS",
                        "APPROVE_STATUS","REMARKS","ROW_NR","COL_NR","WELLTYPE_CODE","COMPOUND_TY","COMPOUND_NR",
                        "CONCENTRATION","IS_VALID","FEATURE_NAME","FEATURE_ALIAS", colname.WELLTYPE_CODE),
  stopIfNoCntrFound = TRUE
)
{
  data.frame.assertColumnsExist(df.welldata, c(groupBy, colname.WELLTYPE_CODE, featureColumns))
  
  cn <- colnames(df.welldata)
  df.norm <- ddply(df.welldata, groupBy, .fun = normalize.generic.singleGroup, colname.WELLTYPE_CODE, WELLTYPE_CODES, featureColumns, func, suffix, colnames.ignore, stopIfNoCntrFound)
  attr(df.norm, 'features.normalized') <- setdiff(colnames(df.norm), cn) # normalized columns = all new columns not in df.welldata before
  return(df.norm)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               normalize (% controls)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
normalize.percentControls <- function(
  df.welldata, 
  groupBy             	= c('BARCODE'),
  colname.WELLTYPE_CODE 	= 'WELLTYPE_CODE',		# name of the column that takes the role of 'well type' (to distinugish controls from samples)
  WELLTYPE_CODES 			= c('LC'),     		 	# WELLTYPE_CODE(S) of the control wells
  featureColumns      	= NULL, 				# columns to normalize, NULL: all numeric columns not in 'colnames.ignore'
  agg.function       	 	= median,
  suffix              	= '.prctCntr',
  colnames.ignore 		= c("EXPERIMENT_ID","EXPERIMENT_NAME","PLATE_ID","BARCODE","SEQUENCE_IN_RUN","PLATE_INFO","VALIDATE_STATUS",
                        "APPROVE_STATUS","REMARKS","ROW_NR","COL_NR","WELLTYPE_CODE","COMPOUND_TY","COMPOUND_NR",
                        "CONCENTRATION","IS_VALID","FEATURE_NAME","FEATURE_ALIAS", colname.WELLTYPE_CODE),
  stopIfNoCntrFound = TRUE
)
{
  return(
    normalize.generic(
      df.welldata, 
      groupBy,
      colname.WELLTYPE_CODE,	# name of the column that takes the role of 'well type' (to distinugish controls from samples)
      WELLTYPE_CODES,      	# WELLTYPE_CODE(S) of the control wells
      featureColumns, 		# columns to normalize, NULL: all numeric columns not in 'colnames.ignore'
      func  			= function(x, cntr) {100 * x / agg.function(cntr)},	# normalize using this function
      suffix,
      colnames.ignore,
      stopIfNoCntrFound
    )
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   z-score normalization with pooled s.d. (substract plate-based mean/median, divide by global and possibly pooled sd (or mad)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
normalize.zscoreControls_plateMean_globalPooledStdev <- function(
  df.welldata, 
  groupBy        			     = NULL,				# these groups will be handled completely independently
  groupBy.mean        	   = c('BARCODE'),		# grouping for means (or medians)
  groupBy.pooled_sd 		   = NULL,   			# grouping for s.d. (or mad), pool the s.d. (or mad) into one value, NULL: s.d. computed across all plates
  colname.WELLTYPE_CODE    = 'WELLTYPE_CODE',
  WELLTYPE_CODES 			     = c('LC'),     	 	# WELLTYPE_CODE(S) of the control wells
  featureColumns      	   = NULL, 			# columns to normalize, NULL: all numeric columns not in 'colnames.ignore'
  suffix              	   = '.zGSCntr',
  robust              	   = FALSE,        	# FALSE: mean / sd, TRUE: median / mad
  centerBeforeComputingSD  = TRUE,				# TRUE: subtract within each "groupBy.mean"-group the mean (or median), then compute s.d. // FALSE: compute s.d. on non-centered data	
  colnames.ignore 		     = c("WELL_ID", "EXPERIMENT_ID","EXPERIMENT_NAME","PLATE_ID","BARCODE","SEQUENCE_IN_RUN","PLATE_INFO","VALIDATE_STATUS",
                        "APPROVE_STATUS","REMARKS","ROW_NR","COL_NR","WELLTYPE_CODE","COMPOUND_TY","COMPOUND_NR",
                        "CONCENTRATION","IS_VALID","FEATURE_NAME","FEATURE_ALIAS")
)
{
  columns_atStart <- colnames(df.welldata)
  nrow_atStart 	<- nrow(df.welldata) 
  
  if(is.null(featureColumns))
  {
    cnames <- colnames(df.welldata)
    cnames <- cnames[!(cnames %in% colnames.ignore)]
  } else {
    cnames <- as.character(featureColumns)
  }
  
  data.frame.assertColumnsExist(df.welldata, c(groupBy, groupBy.mean, groupBy.pooled_sd, colname.WELLTYPE_CODE, cnames))
  
  # get the numeric columns
  cnames <- cnames[as.logical(lapply(df.welldata[, cnames, drop = FALSE], is.numeric))] 
  
  list.df.norm 	<- list()
  cnt 			<- 0
  
  list.df <- dlply(df.welldata, groupBy)
  for(h in 1 : length(list.df)) # loop over independent data groups
  {
    df.welldata <- list.df[[h]]
    
    # subtract mean (or median) per plate (if groupBy.mean is BARCODE)		
    if(centerBeforeComputingSD)
    {
      list.df_h <- dlply(df.welldata, groupBy.mean)
      for(k in 1 : length(list.df_h))  # loop over plates (if group.mean is BARCODE)
      {
        df.welldata_k <- list.df_h[[k]]
        
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        #
        #   get controls wells
        #
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        ix 		<- df.welldata_k[, colname.WELLTYPE_CODE] %in% WELLTYPE_CODES
        x.cntr 	<- as.matrix(df.welldata_k[ix, cnames])
        
        if(nrow(x.cntr) == 0)
        {
          debug.write.csv(df.welldata_k, name = 'No controls')
          stop('No controls found. See Excel table.')
        }
        
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        #
        #   subtract mean (or median)
        #
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        if(robust)
          cen.cntr <- apply(x.cntr, 2, median)
        else
          cen.cntr <- apply(x.cntr, 2, mean)
        
        x <- as.matrix(df.welldata_k[, cnames])
        df.welldata_k[, cnames] <- sweep(x, 2, cen.cntr, '-') 
        
        list.df_h[[k]] <- df.welldata_k 
      }
      df.welldata <- rbind.fill(list.df_h)
    } 
    
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    #     compute pooled s.d. for group h for all features
    #
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    {
      # get all control wells, for global sd/mad
      df.cntr.all <- df.welldata[df.welldata[, colname.WELLTYPE_CODE] %in% WELLTYPE_CODES, ]
      stopifnot(!is.null(df.cntr.all) & nrow(df.cntr.all) >= 6) # at least 6 control wells required for sd / mad. this number is way too low in practice, though! 
      
      cat('groupBy.pooled_sd:', paste0(groupBy.pooled_sd, collapse = ' // '), '\n')
      cat('cnames:', paste0(cnames, collapse = ' // '), '\n')
      
      df.cntrCI <- aggregate.dataframe.3(
        df.cntr.all,                                # the data 
        groupBy     = groupBy.pooled_sd,           	# group data.frame by these variables ==> aggregate each group
        list.agg    = list( 						# list of list(var = myvariables, func = myfunctions, colname = mycolnames), var == NULL ==> use all numeric columns 
          list(                        	
            var     = cnames,    						# aggregate these variables (i.e., columns), NULL: all numeric columns not in c(unique.var, groupBy)
            func    = ifelse(robust, 'mad', 'sd'),      # aggregation functions to be applied (given as strings, not as function handles)    
            colname = c('<agg.var>')    				# naming patterns of the output columns, <agg.var> = placeholder for the resp. input column
          ),
          list(                        	
            var     = cnames[1],    # aggregate these variables (i.e., columns), NULL: all numeric columns not in c(unique.var, groupBy)
            func    = 'length',    	# aggregation functions to be applied (given as strings, not as function handles)    
            colname = 'n__'   		# naming patterns of the output columns, <agg.var> = placeholder for the resp. input column
          )
        )
      )
      
      stopifnot(!is.null(df.cntrCI) & nrow(df.cntrCI) > 0)
      
      if(nrow(df.cntrCI) > 1)
        # compute (n-1) weighted average
      {
        n0 	<- df.cntrCI$n__ - 1
        N 	<- sum(n0) 
        for(i in 1 : length(cnames)) # loop over features
          df.cntrCI[1, cnames[i]] <- sqrt(sum(n0 * df.cntrCI[, cnames[i]]^2) / N)
        
        df.cntrCI <- df.cntrCI[1, ]		# keep only the first row with the pooled s.d. // code below relies on df.cntrCI having only one row!!!	
      }
      df.cntrCI <- df.cntrCI[, !(colnames(df.cntrCI) %in% c('n__', groupBy.pooled_sd)), drop = FALSE] # remove grouping and counts columns
    }
    stopifnot(nrow(df.cntrCI) == 1)
    
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    #   apply mean / sd (or median / mad) on all wells
    #
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    df.welldata <- list.df[[h]]
    list.df_h <- dlply(df.welldata, groupBy.mean)
    for(k in 1 : length(list.df_h))  # loop over plates (if group.mean is BARCODE)
    {
      df.welldata_k <- list.df_h[[k]]
      
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      #
      #   get controls wells
      #
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ix 			<- df.welldata_k[, colname.WELLTYPE_CODE] %in% WELLTYPE_CODES
      df.cntr 	<- df.welldata_k[ix, ]
      
      
      if(nrow(df.cntr) == 0)
      {
        debug.write.csv(df.welldata_k, name = 'No controls')
        stop('No controls found. See Excel table.')
      }
      
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      #
      #  normalize the feature columns (samples and controls)
      #
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      
      for(i in 1:length(cnames)) # loop over features
      {
        colname.i <- cnames[i]
        
        x 		<- 	df.welldata_k[, colname.i]
        cntrCI 	<-  df.cntrCI[, colname.i] # there is only one row in df.cntrCI 
        
        if(robust)
          cntrCI <- 1.4826 * cntrCI # http://en.wikipedia.org/wiki/Median_absolute_deviation
        
        colname.i_norm 	<- paste0(colname.i, suffix)
        
        cntr 			<- df.cntr[,colname.i]
        cntr 			<- cntr[!is.na(cntr) & is.finite(cntr)]
        
        if(length(cntr) == 0)
        {
          debug.write.csv(df.welldata_k, name = 'df_welldata')
          debug.write.csv(df.cntr,  	 name = 'df_cntr')
          stop(sprintf('All controls value are NA or infinite in column %s. See Excel table', colname.i))
        }
        
        if(robust)
          cntr.mid <- median(cntr)
        else
          cntr.mid <- mean(cntr)
        
        x.norm <- (x - cntr.mid) / cntrCI
        
        if(any(is.infinite(x.norm)))
        {
          cat(colname.i,'\n')
          cat('cntr:   ',  paste(cntr, 	collapse = ',\t'),'\n')
          cat('x:      ',  paste(x,     	collapse = ',\t'),'\n')
          cat('x.norm: ',  paste(x.norm, 	collapse = ',\t'),'\n')
          cat('------------------------------------------------\n')
          warning(sprintf('Infinite value produced when normalizing column %s against %s wells. Replacing by NaN\n', colname.i, paste0(WELLTYPE_CODES, collapse = '/')))
          
          x.norm[is.infinite(x.norm)] <- NaN
        }
        df.welldata_k[,colname.i_norm] <- x.norm
      }
      cnt	<- cnt + 1
      list.df.norm[[cnt]] <- df.welldata_k
    }
  }
  
  df.norm <- rbind.fill(list.df.norm)
  attr(df.norm, 'features.normalized') <- setdiff(colnames(df.norm), columns_atStart) # normalized columns = all new columns not in df.welldata before
  
  stopifnot(nrow_atStart == nrow(df.norm))
  return(df.norm)
}


#::::::::::::::::::::::::::::::::::::::::::::::
#
#   pdf tools
#
#::::::::::::::::::::::::::::::::::::::::::::::
init.newPDF <- function(fn.pdf, paper = "a4r", width = 10, height = 7.5, shutDownAllDevices=TRUE, onError.makeNewName = TRUE)
{
  if(shutDownAllDevices)
    dev.off.all()
  
  ensure.dirExists(fn.pdf)
  
  error <- TRUE
  tryCatch ({
    pdf(file = fn.pdf,paper = paper, width = width, height = height)
    error <- FALSE
  },
  error = function(e ) {}
  )
  
  if(error)
  {
    fn.pdf <- file.makeFilenameUnique(fn.pdf, startAt = 1, appendIfAlreadyUnique = FALSE)
    pdf(file = fn.pdf, paper = paper, width = width, height = height)
  }
  return(fn.pdf)
}

dev.off.all <- function()
{
  while(dev.cur() > 1)
    dev.off()
}

finish.pdf <- function(fn.pdf, openPDF = TRUE)
{
  dev.off()
  cat("Graphs written to ",fn.pdf,"\n")
  if(openPDF)
    browseURL(fn.pdf)
}
