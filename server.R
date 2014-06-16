library(shiny)
library(plyr)
library(reshape2)
library(ggplot2)

source("./source/io.R",  local=TRUE)
source("./source/get.R",  local=TRUE)
calculation_dir <- './calculation'
options(shiny.maxRequestSize=30*1024^2)
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
      
      qhts <- input_data_loader()
      if (cytomask) cytoqhts <- cyto_data_loader()
      save_input_curvep(qhts, cytoqhts=cytoqhts, cytomaskthr=cytomaskthr, calculation_dir=calculation_dir)
    })
    
  })
  
  curvep_output_generator <- reactive({
    if(input$run == 0) return(NULL)
    isolate({
      
      #qhts <- input_data_loader()
      #basename <- curvep_input_creator()
      qhts <- curvep_input_creator()
      
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
      
      return(get_curvep_results(qhts, paras, calculation_dir=calculation_dir))
      
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
        qhts <- curvep_output_generator()
        #result <- save_output_data(qhts, logm=logm)
        result <- cbind(qhts$id_hill, qhts$concs, qhts$resps, qhts$curvep_resps)
        write.table(result, file, row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE, append=FALSE)
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
    result <- curvep_output_generator()
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
