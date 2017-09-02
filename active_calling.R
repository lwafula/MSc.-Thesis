source('E:/uhasselt/22yrsem/MThesis/scripts/R_code/stools.R')

# 1. based on the Euclidean distance to the center (= mean) of the DMSO controls each sample is classified as active or inactive
#    - classification threshold is calculated as the 95%-tile of the distribution "distance of DMSO sample to DMSO center"
#    - thus 5% of the DMSO controls are considered as active
# 2. For each treatment, the percentage of active samples (i.e., replicates) is calculated. 
#    (Treatments with >= 50% active replicates are considered active)
#

active_calling <- function(
  datp,
  features,                          # names of the features used in the Euclidean distance, usually the set of selected features from mRMR or another feature selection method  
  active_distance_percentile = 0.95, # classification threshold will be computed as: percentile of the distribution "distance of individual DMSO control to the center of all DMSO controls". 95th-percentile ==> 5% of DMSO controls will be considered active 
  groupBy                    = c('Treatment', 'WELLTYPE_CODE')
)
{
  attribs <- attributes(datp) # restore attributes at the end of the script
  
  #```````````````````````````````````
  #   select columns
  #..................................
  col <- unique(c(groupBy, 'WELL_ID', 'WELLTYPE_CODE', features))
  data <- data.frame.getColumns(datp, col)
  
  
  #```````````````````````````````````````````````````````````````````````````
  #   omit rows containing NA (after selecting required columns!!!)
  #...........................................................................
  nix.NA <- which(!complete.cases(data))
  if (length(nix.NA) > 0)
  {
    warning(paste0(length(nix.NA), ' / ', nrow(data), ' rows containing NA removed.'))
    debug.write.csv(
      data[nix.NA, ],
      name = paste0('removed rows containing NA'),
      na = 'NA'
    )
    data <- data[-nix.NA,]
  }

  #```````````````````````````````````````````````````````````````````````````
  #   define output columns
  #...........................................................................
  cn.distToDMSO 		  <- paste0('Active_distToDMSO')
  cn.active 			    <- paste0('Active_Flag')
  cn.n_active 		    <- paste0('Active_Count')
  cn.prct_activeRepl 	<- paste0('Active_FractionOfRepl')
  cn.n_total 			    <- paste0('ReplicateCount')
  
  #````````````````````````````````````````````
  #   compute activity flag for each sample
  #............................................
  {
    
    #```````````````````````````````````````````
    #    convert to compound-feature-matrix
    #...........................................
    x <- as.matrix(data[, features])
    
    #```````````````````````````````````````````
    #    get control wells
    #...........................................
    nix.DMSO <- which(data$WELLTYPE_CODE == 'LC')
    stopifnot(length(nix.DMSO) > 0)
    
    #```````````````````````````````````````````
    #    compute mean DMSO
    #...........................................
    x_DMSO_mean <- colMeans(x[nix.DMSO, ])
    
    #```````````````````````````````````````````
    #    compute distance to mean DMSO
    #...........................................
    x0 <- sweep(x, 2, x_DMSO_mean, FUN = "-")
    data[, cn.distToDMSO] <- apply(x0, 1, FUN = function(x) sum(x ^ 2)) # sq Euclidean dist to origin (mean DMSO subtracted above)
    
    #```````````````````````````````````````````````````````````````````````````````````````````````````
    #    compute active distance threshold (percentile of "DMSO" replicate - DMSO center distance)
    #...................................................................................................
    thr.distToDMSO <- quantile(data[nix.DMSO, cn.distToDMSO], probs = active_distance_percentile)
    
    #		 hist(log10(data[nix.DMSO, cn.distToDMSO]))
    #		 abline(v = log10(thr.distToDMSO))
    
    #```````````````````````````````````````````
    #    apply threshold
    #...........................................
    data[, cn.active] <- ifelse(data[, cn.distToDMSO] > thr.distToDMSO, 1, 0)
    
    #````````````````````````````````````````````````````````````````````````````
    #   add active flag (and dist2controls) to original welldata
    #............................................................................
    n0 <- nrow(datp)
    c2add <- c(cn.active, cn.distToDMSO)
    datp <-
      datp[, !(colnames(datp) %in% c2add)] # make sure columns to be added are not present in target table
    datp <- merge(datp, data[, c('WELL_ID', c2add)], all.x = TRUE)
    n1 <- nrow(datp)
    stopifnot(n0 == n1)
    stopifnot(!any(is.na(datp$WELL_ID)))
    
    #```````````````````````````````````````````
    #  aggregate active/inactive per treatment
    #...........................................
    list.agg <-
      list(# list of list(var = myvariables, func = myfunctions, colname = mycolnames), var == NULL ==> use all numeric columns
        list(
          # active/inactive
          var     = cn.active,
          # aggregate these variables (i.e., columns), NULL: all numeric columns not in c(unique.var, groupBy)
          func    = c('length', 'sum'),
          # aggregation functions to be applied (given as strings, not as function handles)
          colname = c(cn.n_total, cn.n_active)    				# naming patterns of the output columns, <agg.var> = placeholder for the resp. input column
        ))
    
    df.cmpdAct <- aggregate.dataframe.3(data[, c(groupBy, cn.active)],
                                        groupBy  = groupBy,
                                        list.agg = list.agg)
    
    # fraction of active replicates
    df.cmpdAct[, cn.prct_activeRepl] <- df.cmpdAct[, cn.n_active] / df.cmpdAct[, cn.n_total]
    
  }
  
  stopifnot(!any(is.na(df.cmpdAct$Treatment)))
  tmp <- unique(df.cmpdAct[, groupBy])
  stopifnot(nrow(tmp) == nrow(df.cmpdAct))
  
  #````````````````````````````````````````````````````````````````````````````
  #   add aggr. activity info on treatment level to original welldata
  #............................................................................
  n0 <- nrow(datp)
  c2add <- c(cn.n_total, cn.n_active, cn.prct_activeRepl)
  datp  <- datp[, !(colnames(datp) %in% c2add)] # make sure columns to be added are not present in target table
  
  datp <-
    merge(datp, df.cmpdAct[, c(groupBy, c2add)], by = groupBy, all.x = TRUE)
  n1 <- nrow(datp)
  stopifnot(n0 == n1)
  stopifnot(!any(is.na(datp$WELL_ID)))
  
  #```````````````````````````````````````````
  #   restore attributes of input table
  #...........................................
  datp <- data.frame.setAttributes.noOverwrite(datp, attribs)
  
  #```````````````````````````````````````````
  #   % active treatments
  #...........................................
  n0 <- length(unique(datp$Treatment))
  ix <- datp[, cn.prct_activeRepl] >= 0.5
  n1 <- length(unique(datp$Treatment[ix]))
  cat('% active treatments: ', 100 * (n1 / n0), '\n')
  return(datp)
}

