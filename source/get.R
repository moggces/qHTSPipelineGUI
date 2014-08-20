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
  
  #hill <- qhts[['hills']]
  #identity <- qhts[['id']]
  id_hill <- qhts[['id_hill']]
  concs <- qhts[['concs']]
   
  curvep_result <- read.table(curvep_output, header = TRUE, sep = "\t", quote = "", comment.char = "", skip=skipl, stringsAsFactors = FALSE)
  curvep_paras <- read.table(curvep_output, header = FALSE, sep = "\n", quote = "", comment.char = "", nrow=skipl-2, stringsAsFactors = FALSE)
  
  curvep_resps <- curvep_result[, grepl("Resp[0-9]+", colnames(curvep_result), perl=TRUE)]
  curvep_resps[curvep_resps == -999] <- NA
  colnames(curvep_resps) <- sub('Resp', 'curvep_r', colnames(curvep_resps))
  
  
  id_hill$curvep_wauc <- curvep_result$wAUC
  id_hill$curvep_pod <- curvep_result$POD  
  id_hill$curvep_pod[id_hill$curvep_pod == -999] <- NA
  #hill$curvep_pod = (hill$curvep_pod + 6)*-1 # to be compatible with log10(M)
  
  if (lconc_pod)
  {
    ind <- which(is.na(id_hill$curvep_pod))
    last_concs <- lapply(ind, function (x){
      num <- max(which(! is.na(concs[x, ])))
      last_conc <- concs[x, num]
      return(last_conc)
    })
    id_hill$curvep_pod[ind] <- unlist(last_concs)
  }
  
  id_hill$curvep_remark <- curvep_result$Remarks
  id_hill$curvep_n_corrections <- curvep_result$X.Corrections
  id_hill$curvep_mask <- curvep_result$mask
  id_hill[,"input_mask"] <- qhts$mask
  
  #qhts[['hills']] <- hill
  #qhts[['id_hill']] <- cbind(identity, hill)
  qhts[['id_hill']] <- id_hill
  qhts[['curvep_resps']] <- curvep_resps
  qhts[['paras_lines']] <- curvep_paras # lines not computable
  qhts[['paras']] <- paras 
  qhts[['isinhibitor']] <- isinhibitor
  
  return(qhts)
  
  
}
# to do. average of value < 0 thr  ; value < thr  mask
get_cyto2mask <- function (cytoqhts, thr)
{
  curvep_resps <- cytoqhts$curvep_resps
  if (is.null(curvep_resps)) stop("no curvep resps")
  
  #  if resp is NA, it will become F, F, F, T, NA 
  if (is.na(thr))
  {
    mask <- curvep_resps < 0 
#    apply(curvep_resps, 1, function (x) { })
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

# to include in the front of the output
get_additional_commentl <- function(add_paras)
{
  result <- list()
  if ('cytomaskthr' %in% names(add_paras) ) 
  {
    result[['cytomaskthr']] <- paste('# cytotoxicity response higher than ', 
                                     add_paras[['cytomaskthr']], '% to mask the response.', sep="")
  }
  if ('spiked' %in% names(add_paras) )
  {
    result[['spike_comb']] <- paste('# spike detection (use NCGC outliers).', sep="")
  }
  if ('cytoen' %in% names(add_paras) )
  {
    result[['spike_comb']] <- paste(result[['spike_comb']], '+ cytotoxicity info', sep="")
  }
  if ('uplate' %in% names(add_paras) )
  {
    result[['uplate']] <- paste('# carryover treatment (use plate sequence).', sep="")
  }
  return(result)

}

get_clean_spike <- function (qhts, cytoqhts, lconc_pod)
{
  # collect all the information
  isinhibitor <- qhts$isinhibitor
  thr <- abs(as.numeric(qhts$paras[['thr']]))
 
  # get the responses and mask
  c_respsdf <- qhts$curvep_resps
  respsdf <- qhts$resps
  mask <- qhts$id_hill[, "curvep_mask", drop=FALSE]
  
  # curvep_resp and origianl resp
  if (isinhibitor)
  {
    c_respsdf <- c_respsdf*-1
    respsdf <- respsdf*-1
  }
  
  # cyto curvep resp if available
  if (! is.null(cytoqhts)) cyto_c_respsdf <- cytoqhts$curvep_resps
  
  # the conc with highest resp. if all 0 -> NA , if multiple conc. with same response, the lowest will be reported
  hp <- apply(c_respsdf, 1, function(x) { 
    x <- x[! is.na(x)]
    which(x == max(x) & max(x) != 0  )[1] 
    } )

  spike_ids <- lapply(1:length(hp), function (x) {
      
      spike_id <- vector()
      
      # get the curvep_resps and resps per row
      c_resps <- c_respsdf[x, ! is.na(c_respsdf[x,])]
      resps <- respsdf[x, ! is.na(respsdf[x,])]
      n_concs <- length(c_resps)
      
      # get the cytodata
      cyto_c_resps <- NULL
      if (! is.null(cytoqhts)) cyto_c_resps <- cyto_c_respsdf[x, 1:n_concs] # could be NA
      
      # get the ind of non-outliers
      mask_log <- ! as.logical(as.numeric(strsplit(mask[x,], " ")[[1]]))
      
      # high resp in the middle of the curve
      if (! is.na(hp[x]) &  n_concs - hp[x] > 2 & hp[x] != 1 )
      {
        ind_b <- 1:(hp[x]-1)
        ind_f <- (hp[x]+1):n_concs
        
        ind_b_m <- tail(which(mask_log[ind_b]),1)
        if (length(ind_b_m) == 0) ind_b_m <- hp[x]-1
        ind_f_m <- head(which(mask_log[ind_f]),1)
        if (length(ind_f_m) ==0) ind_f_m <- 1 
        
        # resp before the highest point or response right after the highest point
        hp_max_b <- resps[ind_b_m]
        hp_max_f <- resps[hp[x]+ind_f_m]
        hp_max <- c_resps[hp[x]]
        
        cyto_min_f <- 0
        if ( ! is.null(cyto_c_resps)) cyto_min_f <- min(cyto_c_resps[ind_f], na.rm=TRUE)
        
        #if (hp_max > 150 & hp_max_b < thr & hp_max_f < thr & cyto_min_f == 0) spike_id <- names(hp)[x]
        if ( hp_max_b < thr & hp_max_f < thr & cyto_min_f == 0) spike_id <- names(hp)[x]
        
      }
      return(spike_id)
    }         
  )
  
  # clean the curves
  spike_ids <- unlist(spike_ids)
  if (length(spike_ids) > 0) qhts  <- clean_resps(qhts, spike_ids, lconc_pod, 'spike')
  
  return(qhts)
  
}

get_clean_potent <- function (qhts, lconc_pod)
{
  id_hill <- cbind(qhts$id_hill, qhts$map)
  pot_ids <- which(grepl('INVERSE|POTENT|CARRY_OVER', id_hill$curvep_remark) & id_hill$curvep_wauc != 0)
  
  cols <- c('Library_seq', 'Row', 'Column')
  pot <- id_hill[pot_ids, cols]
  pot <- pot[order(pot$Library_seq),] # from the earlies plate
  
  carry_ids <- vector()
  for (i in 1:nrow(pot))
  {
    l <- pot[i, ] # get the line
    seq <- l$Library_seq
    row <- l$Row
    col <- l$Column
    name <- rownames(l)
    
    if (seq != 1 )
    {
      #pre_id <- id_hill$Library_seq == seq-1 & id_hill$Row == row & id_hill$Column == col
      
      init_seq <- seq
      repeat {
        init_seq <- init_seq - 1
        pre_id <- id_hill$Library_seq ==init_seq  & id_hill$Row == row & id_hill$Column == col
        if (sum(pre_id) != 0 ) break
      } 
      
      pre <- id_hill[pre_id, ] # get the line of previous plate
      pre_act <- pre$curvep_wauc
      pre_name <- rownames(pre)
      pre_mark <- pre$curvep_remark
      first_resp <- qhts$curvep_resps[pre_id, 1]
      if ( (abs(pre_act) > 50 & sum(carry_ids %in% pre_name) == 0) | ( sum(carry_ids %in% pre_name) != 0 & abs(first_resp) > 80 & grepl('INVERSE', pre_mark) ))  
      {carry_ids <- rbind(carry_ids, name) }
    }
  }
  if (length(carry_ids) > 0) qhts  <- clean_resps(qhts, carry_ids, lconc_pod, 'carryover')
  return(qhts)
}

clean_resps <- function (qhts, ids, lconc_pod, comment)
{
  ind <- rownames(qhts$id_hill) %in% ids
  qhts[['id_hill']][ind, 'curvep_remark'] <- paste(qhts$id_hill[ind, 'curvep_remark'], comment, sep="|")
  qhts[['id_hill']][ind, 'curvep_wauc'] <- 0
  qhts[['id_hill']][ind, 'curvep_pod'] <- NA
  if ( lconc_pod) qhts[['id_hill']][ind, 'curvep_pod'] <- apply(qhts$concs[ind, ],1, function (y) { y <- y[!is.na(y)]; tail(y, 1)}) 
  qhts[['curvep_resps']][ind, ] <- 0
  return(qhts)
}
