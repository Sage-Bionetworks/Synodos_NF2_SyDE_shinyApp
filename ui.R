library(shiny)
library(shinyIncubator)

shinyUI( navbarPage("Synodos Data Explorer",
                    
        ##############
        #panel 1
        ##############
        tabPanel("Public Data",
                 progressInit(),
                 h3("Explore NF2 Public Datasets"),
                 sidebarLayout(
                   sidebarPanel(
                     h4('1. Select Grouped Pathways '),
                     selectInput('grouped_pathways', NULL, choices = names(grouped_pathways),
                                 selected=c('NF2_related_pathways'), selectize=T, multiple=T),
                     br(),br(),
                     h4('2. Select individual pathway/s'),
                     selectInput('pathways', NULL, choices = rownames(pathway_enrichment_scores), multiple=TRUE, 
                                 width='400px',selectize=T),
                     br(), br(),
                     h4('3. Select based on genes'), 
                     helpText("Accepts HUGO gene names. Gene names may be separated by comma, space,line, comma "),
                     tags$textarea(id="custom_gene_list",rows=8,cols=400),
                     actionButton("searchByGenes_button", "Search by Genes"),
                     
                     br(),br(), br(),
                    
                     h4("Filtering options "),
                     h5('A. By GEO studies'),
                     selectInput('selected_geo_studies', NULL, choices = unique(publicData_phenotype[['studyId']]),
                                  selectize=T, multiple=T),
                     
                     #condition panel to show available phenotypes
                     conditionalPanel(
                          condition = "output.pubData_nselectedStudies == 1",
                                      br(),
                                      helpText('Available Study Phenotypes'),
                                      selectInput('pubData_selectedPhenotypes', NULL, 
                                                  choices = '',
                                                  selectize=T, multiple=T)
                     ),
                     br(), br(), br(), br(), br(), br(),
                     get_synodos_banner(width="300")
                   ),
                   
                   mainPanel(
                     tabsetPanel(id="pubData_pathway_panel", type="tabs",
                                 tabPanel("Pathway Enrichment",
                                          plotOutput("pubData_pathway_heatmap",height="700px",width="auto",hoverId=NULL),
                                          verbatimTextOutput('pubData_nselectedStudies')
                                 ),
                                 tabPanel("Expression",
                                          plotOutput("pubData_expression_heatmap",height="700px",width="auto",hoverId=NULL)
                                          #verbatimTextOutput('pubData_nselectedStudies')
                                          )
                     )
                   ),
                   fluid=T
                 )
        ),#END TabPanel("Public Data")
        
                  
        # Panel 2: Kinome Screens
        tabPanel("Kinome Screens",
                 h3("Synodos Kinome Screens"),
              sidebarPanel(
                   h4('1. Select Samples'),
                   selectInput('kinome_selected_samples', NULL, choices = c('ALL', as.character(unique(kinomeData$sample))),
                                selected = as.character(unique(kinomeData$sample))[1:2], selectize=T, multiple=T),
                   br(),br(),
                   h4('Filter By:'),
                   br(), br(),
                   
                   h4('2. Select Kinase Family'),
                   selectInput('kinome_selected_kinaseFamily', NULL, choices = c('ALL',unique(kinomeData$Family)),
                               selected = unique(kinomeData$Family)[1:2], selectize=T, multiple=T),
                   br(),br(),br()
                 ),
                 mainPanel(
                   column(3,showOutput("kinome_barPlot", "nvd3"))
                  # column(6,showOutput("kinome_barPlot", "nvd3"))
                   
#                    tabsetPanel(id="kinome_mainPanel", type="tabs",
#                                tabPanel("Explore Data",
#                                         column(10,
#                                                 showOutput("kinome_barPlot", "nvd3")
#                                         )   
#                                ),
#                                tabPanel("Download Data"
#                                       
#                                )
                   
                )
          ),  
        
         # Testing Panel
        tabPanel("testing",
                    sidebarPanel(
                      selectInput(inputId = "x",
                                  label = "Choose X",
                                  choices = c('SepalLength', 'SepalWidth', 'PetalLength', 'PetalWidth'),
                                  selected = "SepalLength"),
                      selectInput(inputId = "y",
                                  label = "Choose Y",
                                  choices = c('SepalLength', 'SepalWidth', 'PetalLength', 'PetalWidth'),
                                  selected = "SepalWidth")
                    ),
                    mainPanel(
                      showOutput("myChart", "polycharts"),
                      showOutput("chart1", "polycharts")
                    )
        ),   
    #navbar pages option      
    fluid = T,
    responsive = T,
    collapsable = T
    )
)

