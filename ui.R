library(shiny)
library(shinyIncubator)

shinyUI( navbarPage("Synodos Data Explorer",
                    
      #panel 1
      #DrugScreen Panel
      tabPanel('Drug Screens',
               h3("Synodos Drug Screens"),
               sidebarPanel(
                
                 
                 h4('1. Select Cell Lines'), 
                 h5('MGH Cell Lines'), 
                 selectInput('MGH_cellLines',NULL, choices = c('ALL', unique(MGH_normViab$cellLine)),
                             selectize=T, multiple=T, selected = c('Syn1', 'Syn2', 'Syn10', 'Syn7')),
                 h5('UCF Cell Lines'),
                 selectInput('UCF_cellLines',NULL, choices = c('ALL',unique(UCF_normViab$cellLine)),
                             selectize=T, multiple=T),
                 tags$a(href=global_cellLines_metadata_link,target="_blank", "cell line metadata"),
                 
                
                 br(),br(),br(),
                 
                 h4('2. Select Drugs'),
                 selectInput('selected_drugs',NULL, choices = c('ALL', unique(drug_normViab$drug)),
                             selectize=T, multiple=T, selected=c('Bortezomib', 'Ganetespib', 'Panobinostat')),
                 tags$a(href=global_drug_metadata_link,target="_blank", "drug metadata"),
                 
                 
                 br(), br(),
                 h4('Plot settings'),
                 selectInput('facet_by','Facet by', choices =c('group', 'experiment', 'cellLine'),
                             selected = c('experiment', 'group'),
                             selectize=T, multiple=T),
                 
                 br(), br()
                 
               ),
               mainPanel(
                 tabsetPanel(id="drug_screens", type="tabs",
                             tabPanel("Viability",
                                      plotOutput("global_drugViab_heatMap",height="700px",width="auto",hoverId=NULL),
                                      br(), br(),
                                      helpText("ps: due to different dynamic range of drug doses across UCF and MGH cellLines
                                               , comparing these together on a heatmap will result in some columns being blank")
                                      
                             ),
                             tabPanel("Max Efficacy",
                                      plotOutput("drug_efficacy",height="700px",width="auto",hoverId=NULL),
                                      helpText("ps: Efficacy is the % of cells that were killed by the given drug")
                             ),
                             tabPanel("ICx",
                                      h4('Select Cell Viability % (ICx)'),
                                      sliderInput('selected_IC_value','IC Value', min=10,
                                                  max=90, step=10, value=50),
                                      plotOutput("drugScreen_ICx_plot",height="700px",width="auto",hoverId=NULL)
                             ),
                             tabPanel("Dose Response",
                                      radioButtons("dose_response_plot_splitBy", label=("split graph by"), 
                                                  choices = c('drug', 'cellLine'), selected = 'drug'),
                                      plotOutput("drugResponse_plots",height="700px",width="auto",hoverId=NULL)
                             )
                 )
               ),
               fluid=T
        ),#END TabPanel("Drug Screen Data")
         
                    
       # Panel 2: Kinome Screens
       tabPanel("Kinome Screens",
                h3("Synodos Kinome Screens"),
                sidebarPanel(
                  h4('1. Select Samples'),
                  selectInput('kinome_selected_samples', NULL, choices=unique(kinomeData$condition),
                              selected = unique(kinomeData$condition)[1], selectize=T, multiple=T),
                  br(),br(),
                  h4('Filter By:'),
                  br(),
                  
                  h4('2. Select Kinase Family'),
                  selectInput('kinome_selected_kinaseFamily', NULL, choices = c('ALL',unique(kinomeData$Family)),
                              selected = unique(kinomeData$Family)[1:2], selectize=T, multiple=T),
                  br(),br(),br()
                ),
                mainPanel(
                  column(3,showOutput("kinome_barPlot", "nvd3"))
                )
       ),  
                    
                    
      #Tab panel 3
      tabPanel("Public Data",
               h3("Explore NF2 Public Datasets"),
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
                        condition = 'input["selected_geo_studies"].length == 1',
                                    br(),
                                    helpText('Available Study Phenotypes'),
                                    selectInput('pubData_selectedPhenotypes', NULL, 
                                                choices = '',
                                                selectize=T, multiple=T)
                   )
                 ),
                 
                 mainPanel(
                   tabsetPanel(id="pubData_pathway_panel", type="tabs",
                               tabPanel("Expression",
                                        plotOutput("pubData_expression_heatmap",height="700px",width="auto",hoverId=NULL)),
                               tabPanel("Pathway Enrichment",
                                        plotOutput("pubData_pathway_heatmap",height="700px",width="auto",hoverId=NULL))
                   )
                 )
      ), #END TabPanel("Public Data")
               
    #navbar pages option      
    fluid = T,
    responsive = T,
    collapsable = T,
    footer = list(helpText('Please report any bugs and/or user experience based feedback at apratap@sagebase.org'),
                  img(src="synodos-banner.jpg", height="150", width="180"))
    )#END nav bar
) #END ShinyUI
