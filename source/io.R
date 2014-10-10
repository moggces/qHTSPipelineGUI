load_input_file <- function (file)
{
  result <- list()
  skipl <- length(grep("^#", readLines(file)))
  
  df <- read.delim(file, quote = "",  skip=skipl, stringsAsFactors=FALSE)
  if (is.null(df$uniqueID)) stop("uniqueID column is required in the user input")
  if (is.null(df$pathway)) stop("pathway column is required in the user input")
  
  #rownames(df) <- df$uniqueID
  rownames(df) <- paste('u', 1:nrow(df))
  
  ### SORT!!! unique ID ##### don't change the sequence order 
  df <- df[order(df$pathway, df$uniqueID),]
  
  
  if (! is.null(df$parent)) df$parent[is.na(df$parent)] <- ''
  
  col_names <- colnames(df)
  
  # the paramters from either Curve Class 
  hill_cols <- c('Curve.Class','Curve.Class2', 'LogAC50', 'Hill.Coef', 'Inf.Activity', 'Zero.Activity', 'Mask.Flags', 'R2', 'p.value')
  hill_cols <- intersect(col_names, hill_cols)
  
  map_cols <- c('Row', 'Column')
  map_cols <- intersect(col_names, map_cols)
  
  x_cols <- grep("conc[0-9]+", col_names, value = TRUE)
  y_cols <- grep("resp[0-9]+", col_names, value = TRUE)
  curvepr_cols <- grep("curvep_r[0-9]+", col_names, value = TRUE)
  
  # mask column
  if (is.null(df$mask) ) df[, "mask"] <- ''
  if (sum(is.na(df$mask)) == nrow(df)) df[, "mask"] <- ''
  
  result[['concs']] <- df[, x_cols, drop=FALSE]
  result[['resps']] <- df[, y_cols, drop=FALSE]
  result[['hills']] <- df[, hill_cols, drop=FALSE]
  result[['curvep_resps']] <- df[, curvepr_cols, drop=FALSE]
  result[['map']] <- df[, map_cols, drop=FALSE]
  result[['filen']] <- file
  result[['mask']] <- df[, "mask", drop=FALSE]
  result[['id']] <- df[, ! colnames(df) %in% c(x_cols, y_cols, hill_cols, map_cols, curvepr_cols, "mask"), drop=FALSE]
  result[['id_hill']] <- cbind(result$id, result$hills)
  return(result)
  
}


save_input_curvep <- function (qhts, cytoqhts, cytomaskthr=NA, spiked=FALSE, calculation_dir)
{
  basename <- as.numeric(as.POSIXct(Sys.time()))
  basename <- as.character(basename)
  basename <- strsplit(basename, ".", fixed=TRUE)[[1]][1]
  
  output <- paste(calculation_dir, "/",  basename, ".htsx", sep="")
  output_p <- paste(output, sep="")
  
  if (cytomaskthr == '') cytomaskthr <- NA
  cytomaskthr <- as.numeric(cytomaskthr)
  
  identity <- qhts[['id']]
  resps <- qhts[['resps']]
  concs <- qhts[['concs']]
  hills <- qhts[['hills']]
  maskdf <- qhts[['mask']]
  
  maskdf$uniqueID <- identity$uniqueID
  
  if (! is.null(cytoqhts)) 
  {
    cytoqhts <- get_cyto2mask(cytoqhts, cytomaskthr)
    cytomaskdf <- cytoqhts$mask
    cytomaskdf[, "uniqueID"] <- cytoqhts[['id']]$uniqueID
    maskdf <- join(subset(maskdf, select=-mask), cytomaskdf, by="uniqueID")
  }
  
  # for ease, Mask.Flags is used when spiked is used and mask column is all blank ''
  if (spiked)
  {
    if ( sum(maskdf$mask == '') == nrow(maskdf)) maskdf[, 'mask'] <- hills$Mask.Flags
  }

  mask <- maskdf$mask
  qhts[['mask']] <- subset(maskdf, select=mask)

  resps[is.na(resps)] <- -999
  concs[is.na(concs)] <- -999 

  #id <- seq(nrow(concs))
  id <- identity$uniqueID
  
  write.table(cbind(nrow(resps), ncol(resps)), output_p,  row.names = FALSE, col.names = FALSE, sep="\t", quote=FALSE, append=FALSE)
  suppressWarnings(write.table(cbind(id, mask, concs, resps), output_p,  row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE, append=TRUE))
  
  if ( ! file.exists(output_p) ) basename <- NULL
  qhts[['basename']] <- basename
  return (qhts)
}

# save_output_data <- function(qhts, logm)
# {
#   concs <- qhts$concs
#   if (logm) concs <- log10(concs)
#   result <- cbind(qhts$id, qhts$hills, qhts$concs, qhts$resps, qhts$curvep_resps)
#   return(result)
# }