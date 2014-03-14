
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Curvep interface"),
  sidebarPanel(
    
    h4('Data'),
    fileInput('file1', 'Import a file with required columns', multiple=FALSE),
    helpText("Required columns: uniqueID, conc[0-9], resp[0-9] (-100% ~ 100%, 0%: baseline). Optional: Mask.Flags"),
    
    tags$br(),
    
    fileInput('cytofile', 'Import a cytotoxicity file with required columns', multiple=FALSE),
    helpText("Required columns: curvep_r[0-9] if mask generation"),
    selectInput("cunit", "concentration unit", 
                choices = list('lgM'='lgM', 'nM' = 'nM', 'uM'='uM', 'M'='M')),
  
    tags$hr(),
    
    h4('Primary parameters'),
    textInput('mxdv', 'max.deviation: ', '5'),
    textInput('thr', 'baseline threshold: ', '15'),
    textInput('rng', 'max.response: ', '-100'),
    checkboxInput("lconc_pod", "use last concentration as default for POD or NA", TRUE),
  
    tags$br(),
    
    h4('Secondary parameters'),
    textInput('cro', 'carryover threshold: ', '80'),
    textInput('ushape', 'min.#points for u-shape: ', '2'),
    textInput('bshift', 'min.#flat points to detect baseline-shift: ', '3'),
    checkboxInput("cort_bshift", "baseline shift correction mode", FALSE),
    checkboxInput("bylo", "favors corrections based on low conc-s", TRUE),
    #checkboxInput("xtinf", "outputs additional curve metrics", TRUE),
    checkboxInput("xplax", "allow extrapolation beyond test conc. boundaries", FALSE),
    checkboxInput("cytomask", "use cytotoxicity data to generate mask", FALSE),
    
    tags$br(),
    
    h4('Additional curve curation'),
    tags$br(),
    
    h5('Spike detection'),
    checkboxInput("ucyto", "use cytotoxicity data", FALSE),
    checkboxInput("uloose", "use loose", FALSE),
    
    tags$br(),
    
    h5('Carryover detection'),
    checkboxInput("uplate", "use plate sequence", FALSE),
    
    tags$br(),
    
    actionButton('run', 'Run'),
    
    tags$hr(),
    downloadButton('save', 'Save')
    #checkboxInput("logM", "concentration is log10(M)", TRUE)
  ),  
 
    
  mainPanel(
    
    tabsetPanel(
      tabPanel( 'Test1', dataTableOutput('test1')),
      tabPanel( 'Test2', textOutput('test2'))
      #tabPanel( "Plot", plotOutput("plot", height="auto", width="auto"))
      #tabPanel("Data", dataTableOutput('qhts_data')),
      #tabPanelAbout()
    )
  )
    
    
  
  
  
))