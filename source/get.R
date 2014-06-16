run_curvep <- function (basename, paras, calculation_dir)
{
  #curvep_exe <- './calculation/curvep'
  curvep_exe <- paste(calculation_dir, "/curvep", sep="")
  curvep_input <- paste( calculation_dir, "/", basename, '.htsx', sep="")
  f_cor <- '-BYLO'
  #extra <- ''
  cort_bshift_v <- ''
  xplax_c <- ''
  
  if (! paras$bylo) f_cor <- '-BYHI'
  #if (paras$xtinf) extra <- '-XTINF'
  if (paras$cort_bshift) cort_bshift_v <- '-BLFX'
  if (paras$xplax) xplax_c <- '-XPLAX'
  
  err <- system(paste(curvep_exe, " ", curvep_input, 
                      " -RNG=", paras$rng, 
                      " -THR=", paras$thr, 
                      " -MXDV=", paras$mxdv, 
                      " -CRO=", paras$cro,
                      " -HTSX=", paras$cunit,
                      " -USHAPE=", paras$ushape, 
                      " -BSHIFT=", paras$bshift,
                      " ", f_cor, 
                      " ", cort_bshift_v,
                      " ", xplax_c, " -LOG",
                      sep=""), 
                intern=TRUE, ignore.stdout = FALSE, ignore.stderr=FALSE)
  return(err)

}

get_curvep_results <- function (qhts,  paras, calculation_dir) 
{ 
  basename <- qhts$basename
  
  err <- run_curvep(basename, paras, calculation_dir)
  if (err !=  "Done.") stop (err)
  lconc_pod <- paras$lconc_pod
  curvep_output <- paste( calculation_dir, "/", basename, "-mdf.htsx", sep="")
  isinhibitor <- TRUE
  if (paras$rng > 0 ) isinhibitor <- FALSE
  
  
  skipl <- length(grep("^#", readLines(curvep_output))) + 1 # the number of row number of conc
  hill <- qhts[['hills']]
  identity <- qhts[['id']]
  concs <- qhts[['concs']]
   
  curvep_result <- read.table(curvep_output, header = TRUE, sep = "\t", quote = "", comment.char = "", skip=skipl)
  curvep_resps <- curvep_result[, grepl("Resp[0-9]+", colnames(curvep_result), perl=TRUE)]
  curvep_resps[curvep_resps == -999] <- NA
  colnames(curvep_resps) <- sub('Resp', 'curvep_r', colnames(curvep_resps))
  
  hill$curvep_wauc <- curvep_result$wAUC
  hill$curvep_pod <- curvep_result$POD  
  hill$curvep_pod[hill$curvep_pod == -999] <- NA
  #hill$curvep_pod = (hill$curvep_pod + 6)*-1 # to be compatible with log10(M)
  
  if (lconc_pod)
  {
    ind <- which(is.na(hill$curvep_pod))
    last_concs <- lapply(ind, function (x){
      num <- max(which(! is.na(concs[x, ])))
      last_conc <- concs[x, num]
      return(last_conc)
    })
    hill$curvep_pod[ind] <- unlist(last_concs)
  }
  
  hill$curvep_remark <- curvep_result$Remarks
  hill$curvep_n_corrections <- curvep_result$X.Corrections
  hill$curvep_mask <- curvep_result$mask
  hill[,"input_mask"] <- qhts$mask
  
  qhts[['hills']] <- hill
  qhts[['id_hill']] <- cbind(identity, hill)
  qhts[['curvep_resps']] <- curvep_resps
  qhts[['paras']] <- paras
  qhts[['isinhibitor']] <- isinhibitor
  
  return(qhts)
  
  
}


### get the valid responses after curvep. the next concentration of the valid responses
# get_cyto2mask <- function (cytoqhts)
# {
#   curvep_resps <- cytoqhts$curvep_resps
#   if (is.null(curvep_resps)) stop("no curvep resps")
#   mask <- curvep_resps < 0
#   max_n_conc <- ncol(cytoqhts[['concs']])
#   max.n <- vector()
#   
#   mask <- apply(mask, 2, as.numeric)
#   l <- apply(mask, 1, function(x) { if ( sum(x == 1, na.rm=TRUE)  > 0 ) {  min(which(x == 1)) } else { NA } })
#   
#   for ( x in 1:nrow(mask) )
#   {
#     if ( ! is.na(l[x]) )
#     {
#       if (l[x] == 1 )
#       {
#         mask[x, ] <- 0
#       } else if ( l[x] != max_n_conc )
#       {
#         mask[x, l[x]] <- 0
#       } else if ( l[x] == max_n_conc )
#       {
#         mask[x, ] <- 0
#       }
#     }
#   }
#   mask.n <- apply(mask, 1, paste, collapse=" ")
#   return (mask.n)
# }



get_cyto2mask <- function (cytoqhts, thr)
{
  curvep_resps <- cytoqhts$curvep_resps
  if (is.null(curvep_resps)) stop("no curvep resps")
  
  if (is.na(thr))
  {
    mask <- curvep_resps < 0 
  } else
  {
    mask <- curvep_resps < thr*-1 ##
  }
  
  #max_n_conc <- ncol(cytoqhts[['concs']])
  #max.n <- vector()
  
  mask <- apply(mask, 2, as.numeric) # columnwise
  
  # find the lowest conc that is smaller than threshold
  l <- apply(mask, 1, function(x) { if ( sum(x == 1, na.rm=TRUE)  > 0 ) {  min(which(x == 1)) } else { NA } })
  
  mask <- do.call("rbind", lapply(1:nrow(mask), 
        function (x)
        {
          if (! is.na(l[x]))
          {
            if (l[x] == 1) mask[x, ] <- 0
            if (is.na(thr))
            {
              mask[x, l[x]] <- 0
            }
          }
          return(mask[x,])
        }))
  

  mask <- apply(mask, 1, paste, collapse=" ")
  cytoqhts[["mask"]][, "mask"] <- mask
  return (cytoqhts)
}