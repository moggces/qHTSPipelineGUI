load_input_file <- function (file)
{
  result <- list()
  
  df <- read.delim(file, quote = "", stringsAsFactors=FALSE)
  if (is.null(df$uniqueID)) stop("uniqueID column is required in the user input")
  if (is.null(df$pathway)) stop("pathway column is required in the user input")
  df <- df[order(df$pathway, df$uniqueID),]   ### SORT!!! unique ID
  if (! is.null(df$parent)) df$parent[is.na(df$parent)] <- ''
  
  col_names <- colnames(df)
  
  hill_cols <- c('Curve.Class','Curve.Class2', 'LogAC50', 'Hill.Coef', 'Inf.Activity', 'Zero.Activity', 'Mask.Flags', 'R2', 'p.value')
  hill_cols <- intersect(col_names, hill_cols)
  
  map_cols <- c('Row', 'Column')
  map_cols <- intersect(col_names, map_cols)
  
  x_cols <- grep("conc[0-9]+", col_names, value = TRUE)
  y_cols <- grep("resp[0-9]+", col_names, value = TRUE)
  curvepr_cols <- grep("curvep_r[0-9]+", col_names, value = TRUE)
  
  result[['concs']] <- df[, x_cols, drop=FALSE]
  result[['resps']] <- df[, y_cols, drop=FALSE]
  result[['hills']] <- df[, hill_cols, drop=FALSE]
  result[['curvep_resps']] <- df[, curvepr_cols, drop=FALSE]
  result[['map']] <- df[, map_cols, drop=FALSE]
  result[['filen']] <- file
  result[['id']] <- df[, ! colnames(df) %in% c(x_cols, y_cols, hill_cols, map_cols, curvepr_cols), drop=FALSE]
  
  return(result)
  
}


save_input_curvep <- function (qhts, cytoqhts, calculation_dir)
{
  basename <- as.numeric(as.POSIXct(Sys.time()))
  basename <- as.character(basename)
  basename <- strsplit(basename, ".", fixed=TRUE)[[1]][1]
  
  output <- paste(calculation_dir, "/",  basename, ".htsx", sep="")
  output_p <- paste(output, sep="")
  
  identity <- qhts[['id']]
  resps <- qhts[['resps']]
  concs <- qhts[['concs']]
  hills <- qhts[['hills']]
  
  hills$pathway <- identity$pathway
  
  if (is.null(hills$Mask.Flags)) hills$Mask.Flags <- ''
  
  if (! is.null(cytoqhts)) 
  {
    mask <- get_cyto2mask(cytoqhts)
    p <- hills$pathway
    l <- split(hills, p)
    l <- lapply(l, transform, Mask.Flags = mask )
    hills <- unsplit(l, p)
  }

  mask <- hills$Mask.Flags
  
  resps[is.na(resps)] <- -999
  
#   if (toUM)
#   {
#     concs <- 10^concs*1000000  # has to be not log value and uM in current version
#     concs[is.na(concs)] <- -999
#   }
  concs[is.na(concs)] <- -999 

  id <- seq(nrow(concs))
  
  write.table(cbind(nrow(resps), ncol(resps)), output_p,  row.names = FALSE, col.names = FALSE, sep="\t", quote=FALSE, append=FALSE)
  suppressWarnings(write.table(cbind(id, mask, concs, resps), output_p,  row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE, append=TRUE))
  
  if ( ! file.exists(output_p) ) basename <- NULL
    
  return (basename)
}

# save_output_data <- function(qhts, logm)
# {
#   concs <- qhts$concs
#   if (logm) concs <- log10(concs)
#   result <- cbind(qhts$id, qhts$hills, qhts$concs, qhts$resps, qhts$curvep_resps)
#   return(result)
# }