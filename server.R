library(shiny)
library(plyr)
library(reshape2)
library(ggplot2)

source("./source/io.R",  local=TRUE)
source("./source/get.R",  local=TRUE)
calculation_dir <- './calculation'
options(shiny.maxRequestSize=30*1024^2)
#known bugs: if the cytomask generation. if NA , it will have NA in the mask
############
shinyServer(function(input, output) {
   
  input_data_loader <- reactive({
    inFile <- input$file1
    
    if(input$run == 0) return(NULL)
    
    if (is.null(inFile)) return(NULL)
    
    isolate({
      path <- inFile$datapath
      load_input_file(path)
    }) 
  })
  
  cyto_data_loader <- reactive({
    inFile <- input$cytofile
    
    if(input$run == 0) return(NULL)
    
    if (is.null(inFile)) return(NULL)
    
    isolate({
      path <- inFile$datapath
      load_input_file(path)
    }) 
  })
  
  curvep_input_creator <- reactive({
    if(input$run == 0) return(NULL)
    
    isolate({
      cytomask <- input$cytomask
      cytoqhts <- NULL
      cytomaskthr <- input$cytomaskthr
      spiked <- input$spiked
      
      if (cytomask & spiked) stop("only one type of mask is allowed")
      
      qhts <- input_data_loader()
      if (cytomask) cytoqhts <- cyto_data_loader()
      
      save_input_curvep(qhts, cytoqhts=cytoqhts, cytomaskthr=cytomaskthr, spiked=spiked, calculation_dir=calculation_dir)
    })
    
  })
  
  curvep_output_generator <- reactive({
    if(input$run == 0) return(NULL)
    isolate({
      
      #qhts <- input_data_loader()
      #basename <- curvep_input_creator()
      qhts <- curvep_input_creator()
      
      # r specific parameters
      add_paras <- list()
      more_comments <- list()
      
      cytomask <- input$cytomask
      cytomaskthr <- input$cytomaskthr
      spiked <- input$spiked
      cytoen <- input$cytoen
      uplate <- input$uplate
      
      if (cytomask) add_paras[['cytomaskthr']] <- cytomaskthr
      if (uplate) add_paras[['uplate']] <- uplate
      if (spiked) 
      {
        add_paras[['spiked']] <- spiked
        if (cytoen) add_paras[['cytoen']] <- cytoen
      }
      if ( length(add_paras) > 0 ) more_comments <- get_additional_commentl(add_paras)
     
      # curvep specific parameters
      mxdv <- input$mxdv
      thr <- input$thr
      rng <- input$rng
      cro <- input$cro
      cort_bshift <- input$cort_bshift
      ushape <- input$ushape
      bshift <- input$bshift
      bylo <- input$bylo
      #xtinf <- input$xtinf
      lconc_pod <- input$lconc_pod
      cunit <- input$cunit
      xplax <- input$xplax
      
      paras <- list(rng=rng, thr=thr, mxdv=mxdv, cro=cro, ushape=ushape, bshift=bshift, 
                    cort_bshift=cort_bshift, bylo=bylo, lconc_pod=lconc_pod, cunit=cunit, xplax=xplax)
      
      result <- get_curvep_results(qhts, paras, calculation_dir=calculation_dir)
      result[['paras_lines']] <- as.data.frame(unlist(c(result[['paras_lines']], more_comments)))
      return(result)
      
    })
  })
  
  post_curvep_curation <- reactive({
    if(input$run == 0) return(NULL)
    isolate({
      spiked <- input$spiked
      cytoen <- input$cytoen
      lconc_pod <- input$lconc_pod
      uplate <- input$uplate
      cytoqhts <- NULL
      qhts <- curvep_output_generator()
      if (cytoen & spiked) cytoqhts <- cyto_data_loader()
      if (spiked) qhts <- get_clean_spike(qhts, cytoqhts, lconc_pod)
      if (uplate) qhts <- get_clean_potent(qhts, lconc_pod)
      return (qhts)
      
    })
  })
  
  output$save <-  downloadHandler(
    filename = function() {
      inFile <- input$file1
      if (is.null(inFile)) return(NULL)
      basename <- sub("(.*)\\..*", "\\1", inFile$name) #sub("(.*\\/)([^.]+)(\\.[[:alnum:]]+$)", "\\2", path)
      filename <- paste(basename, "_curvep.txt", sep="")
      return(filename)
    },
    content = function(file) {
      if(input$run == 0) return(NULL)
      isolate({
        #logm <- input$logM
        #qhts <- curvep_output_generator()
        qhts <- post_curvep_curation()
        #result <- save_output_data(qhts, logm=logm)
        result <- cbind(qhts$id_hill, qhts$concs, qhts$resps, qhts$curvep_resps,qhts$map)
        
        # also add the comment lines
        write.table(qhts$paras_lines, file, row.names = FALSE, col.names = FALSE, sep="\t", quote=FALSE,append=FALSE)
        write.table(result, file, row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE, append=TRUE)
      })
    }
  )
  
  output$test1 <- renderDataTable({
    if(input$run == 0) return(NULL)
    #input_data_loader()$id
    #as.data.frame(curvep_input_creator())
    #bb <- curvep_input_creator()
    #strsplit(bb, ".", fixed=TRUE)[[1]][1]
    #bb
    result <- post_curvep_curation()
    return(cbind(result$id_hill, result$resps, result$curvep_resps))
  })
  
  output$test2 <- renderText({
    if(input$run == 0) return(NULL)
    
    return(as.character(input$thr))
    #input_data_loader()$id
    #as.data.frame(curvep_input_creator())
    #bb <- curvep_input_creator()
    #strsplit(bb, ".", fixed=TRUE)[[1]][1]
    #bb
    #result <- curvep_output_generator()
    #result[['concs']]
  })
  
})
