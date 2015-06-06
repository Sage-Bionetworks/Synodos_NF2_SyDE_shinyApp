library(shiny)

shinyUI( navbarPage("Synodos Data Explorer",    
  #DrugScreen Panel 1
  tabPanel('Drug Screens',
     bsCollapse(multiple = F, id = "drugScreen", open="Collapse",
        bsCollapsePanel("Collapse",
           fluidRow(
              column(width=3, 
                 h5('1. Select MGH Cell Lines'),
                 selectInput('MGH_cellLines',NULL, choices = c('ALL', unique(MGH_normViab$cellLine)),
                             selectize=T, multiple=T, selected = c('Syn1', 'Syn2', 'Syn10', 'Syn7')),
                 tags$a(href=global_cellLines_metadata_link,target="_blank", "cell line metadata")),
              column(width=3,
                 h5('2. Select UCF Cell Lines'),
                 selectInput('UCF_cellLines',NULL, choices = c('ALL',unique(UCF_normViab$cellLine)),
                             selectize=T, multiple=T)),
              column(width=3,
                 h5('3. Select Drugs'),
                 selectInput('selected_drugs',NULL, choices = c('ALL', unique(drug_normViab$drug)),
                              selectize=T, multiple=T, selected=c('Bortezomib', 'Ganetespib', 'Panobinostat')),
                 tags$a(href=global_drug_metadata_link,target="_blank", "drug metadata")),
              column(width=3,
                     h5('4. Plot Settings'),
                     selectInput('facet_by','Facet by', choices =c('group', 'experiment', 'cellLine'),
                                  selected = c('experiment', 'group'),
                                  selectize=T, multiple=T),
                     submitButton("Update") ) 
           ) #End FluidRow
        ) #End BS Collapse Panel
    ), # end BS Collapse
    hr(),
    tabsetPanel(id="drug_screens", type="tabs",
      tabPanel("Viability",
          checkboxInput("drugViability_heatmap_col_cluster", label=("cluster drugs (y-axis)"), value=TRUE),
          plotOutput("global_drugViab_heatMap",height="600px",width="auto",hoverId=NULL),
          br(), br(),
          helpText("The heatmap above compares the normalized cell viablity(%) across the drug dosages(x-axis) \n
                      Each y-axis row is a replicate for the selected cellLines and drugs."),                                    
          helpText("ps: due to different dynamic range of drug doses across UCF and MGH cellLines
                    ,comparing these together on a heatmap will result in some columns being blank")),                             
          tabPanel("Max Efficacy",
                    plotOutput("drug_efficacy",height="600px",width="auto",hoverId=NULL),
                    helpText("ps: Efficacy is the % of cells inhibhited by the drug")),
          tabPanel("ICx",
                  h4('Select Cell Viability % (ICx)'),
                  sliderInput('selected_IC_value','IC Value', min=10, max=90, step=10, value=50),
                  plotOutput("drugScreen_ICx_plot",height="700px",width="auto",hoverId=NULL)),
          tabPanel("Dose Response",
                  radioButtons("dose_response_plot_splitBy", label=("split graph by"), 
                              choices = c('drug', 'cellLine'), selected = 'drug'),
                  plotOutput("drugResponse_plots",height="700px",width="auto",hoverId=NULL))
    ) #END tabset panel
  ), #END drug screen UI tab panel

###############################################                    

# Panel 2: Kinome Screens
 tabPanel("Kinome Screens", open = T,
  bsCollapse(multiple = T, id = "kinome", open="Select Samples",
   bsCollapsePanel("Select Samples", 
    fluidRow( column(width=3, 
                h5('1. Select Cell Line'),
                selectInput('kinome_selected_cellLines', NULL, choices=unique(kinomeData$cellLine),
                             selected = unique(kinomeData$cellLine)[1], selectize=T, multiple=T),
                tags$a(href=global_cellLines_metadata_link,target="_blank", "cell line metadata")
              ), #end column  
              column(width=3,
                h5('2. Drug'),
                selectInput('kinome_selected_drugs', NULL, choices=unique(kinomeData$drug),
                            selected = unique(kinomeData$drug)[1], selectize=T, multiple=T)
              ), #end column
              column(width=3,
                h5('3. replicate'),
                selectInput('kinome_selected_replicates', NULL, choices=unique(kinomeData$replicate),
                            selected = unique(kinomeData$replicate)[1], selectize=T, multiple=T)
              ), #end column
              column(width=3,
                h5('4. time'),
                selectInput('kinome_selected_time', NULL, choices=unique(kinomeData$time),
                            selected = unique(kinomeData$time)[1], selectize=T, multiple=T),
                submitButton("Update") 
              ) #end column
    )# end fluid row
  ),  # end bascollapse panel
  bsCollapsePanel("Filters",
    fluidRow(
      column(width=3, h5('1. Select Kinase Family'),
             selectInput('kinome_selected_kinaseFamily', NULL, choices = c('ALL',unique(kinomeData$Family)),
                         selected = c('AGC', 'CAMK'), selectize=T, multiple=T)),
      column(width=3, h5('2. Select Genes'),
             selectInput('kinome_selected_genes', NULL, choices = c(unique(kinomeData$Gene)),
                         selectize=T, multiple=T)),
      column(width=3,sliderInput("kinome_var_threshold", h5("3. % variation"), min=0, max=200,   
                           value=c(0,200), step=5)),
      column(width=3,sliderInput('kinome_ratio_includeRegion', h5('4. Include Ratio between:'),
                                  min=1, max=round(max(abs(kinomeData$ratio))), value=c(1,round(max(abs(kinomeData$ratio)))), 
                                  step=0.1))
    ), # end fluid row
    fluidRow(
      column(width=3,offset=6,sliderInput("kinome_PeptideCount_threshold", h5("5. Peptide Counts"), 
                                  min=0, max=max(kinomeData$count), value=c(1,max(kinomeData$count)),step=1)),
      column(width=3,submitButton("Update") )
    ) # end fluid row               
  ) # end bsCollapse panel - Filters
 ),  # end bscollapse
  hr(),
 #Kinome Display
  fluidRow(column(width=6, h5('Protein Ratio(log2) barplot'),
                  plotOutput("kinome_barPlot")),
           column(width=6, br(), h5('Scatterplot'), 
                  rbokehOutput('kinome_scatterPlot'))),
           
   fluidRow(column(width=4, offset=6, 
                   br(),
                   h5('Kinome Ratio Histogram'),
                   plotOutput('kinome_ratio_histogram')))
 ), #end kinome tab panel


#       #Tab panel 3
#       tabPanel("Public Data",
#         sidebarPanel(
#                     h4('1. Select Grouped Pathways '),
#                     selectInput('grouped_pathways', NULL, choices = names(grouped_pathways),
#                                selected=c('NF2_related_pathways'), selectize=T, multiple=T),
#                     br(),br(),
#                     h4('2. Select individual pathway/s'),
#                     selectInput('pathways', NULL, choices = rownames(pathway_enrichment_scores), multiple=TRUE, 
#                                width='400px',selectize=T),
#                     br(), br(),
#                     h4('3. Select based on genes'), 
#                     helpText("Accepts HUGO gene names. Gene names may be separated by comma, space,line, comma "),
#                     tags$textarea(id="custom_gene_list",rows=8,style="width:300px"),
#                     actionButton("searchByGenes_button", "Search by Genes"),
#                    
#                     br(),br(), br(),
#                 
#                     h4("Filtering options "),
#                     h5('A. By GEO studies'),
#                     selectInput('selected_geo_studies', NULL, choices = unique(publicData_phenotype[['studyId']]),
#                                 selectize=T, multiple=T),
#                  
#                    #condition panel to show available phenotypes
#                    conditionalPanel(
#                         condition = 'input["selected_geo_studies"].length == 1',
#                                     br(),
#                                     helpText('Available Study Phenotypes'),
#                                     selectInput('pubData_selectedPhenotypes', NULL, 
#                                                 choices = '',
#                                                 selectize=T, multiple=T)
#                  )
#                ),
#                mainPanel(
#                  tabsetPanel(id="pubData_pathway_panel", type="tabs",
#                              tabPanel("Expression",
#                                       plotOutput("pubData_expression_heatmap",height="700px",width="auto",hoverId=NULL)),
#                              tabPanel("Pathway Enrichment",
#                                       plotOutput("pubData_pathway_heatmap",height="700px",width="auto",hoverId=NULL))
#                  )
#                )     
#       ), #END TabPanel("Public Data")
               
#navbar pages option      
 fluid = T,
 responsive = T,
 collapsible = T
#  footer = list(img(src="synodos-banner.jpg", height="60", width="120"),
#                 helpText('Please report any bugs and/or user experience based feedback at apratap@sagebase.org')
#                )
  )#END nav bar
) #END ShinyUI
