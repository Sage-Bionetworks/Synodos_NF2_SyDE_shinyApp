library(shiny)
library(shinyIncubator)
library(shinyBS)



shinyUI( navbarPage("Synodos Data Explorer",
                    
      #panel 1
      #DrugScreen Panel
      tabPanel('Drug Screens',
               bsCollapse(multiple = F, id = "drugScreen", open="drugScreen1",
                          bsCollapsePanel("Collapse", id='drugScreen1',
                             fluidRow(
                                column(width=3, 
                                       h5('1. Select MGH Cell Lines'),
                                       selectInput('MGH_cellLines',NULL, choices = c('ALL', unique(MGH_normViab$cellLine)),
                                                   selectize=T, multiple=T, selected = c('Syn1', 'Syn2', 'Syn10', 'Syn7')),
                                       br(),
                                       tags$a(href=global_cellLines_metadata_link,target="_blank", "cell line metadata")
                                ),
                                column(width=3,
                                       h5('2. Select UCF Cell Lines'),
                                       selectInput('UCF_cellLines',NULL, choices = c('ALL',unique(UCF_normViab$cellLine)),
                                                   selectize=T, multiple=T)           
                                 ),
                                column(width=3,
                                       h5('3. Select Drugs'),
                                       selectInput('selected_drugs',NULL, choices = c('ALL', unique(drug_normViab$drug)),
                                                    selectize=T, multiple=T, selected=c('Bortezomib', 'Ganetespib', 'Panobinostat')),
                                       tags$a(href=global_drug_metadata_link,target="_blank", "drug metadata")
                                ),
                                column(width=3,
                                       h5('4. Plot Settings'),
                                       selectInput('facet_by','Facet by', choices =c('group', 'experiment', 'cellLine'),
                                                    selected = c('experiment', 'group'),
                                                    selectize=T, multiple=T)
                                ) 
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
                                              ,comparing these together on a heatmap will result in some columns being blank")      
                          ),                             
                          tabPanel("Max Efficacy",
                                    plotOutput("drug_efficacy",height="600px",width="auto",hoverId=NULL),
                                    helpText("ps: Efficacy is the % of cells inhibhited by the drug")
                          ),
                          tabPanel("ICx",
                                  h4('Select Cell Viability % (ICx)'),
                                  sliderInput('selected_IC_value','IC Value', min=10, max=90, step=10, value=50),
                                  plotOutput("drugScreen_ICx_plot",height="700px",width="auto",hoverId=NULL)
                          ),
                          tabPanel("Dose Response",
                                  radioButtons("dose_response_plot_splitBy", label=("split graph by"), 
                                              choices = c('drug', 'cellLine'), selected = 'drug'),
                                  plotOutput("drugResponse_plots",height="700px",width="auto",hoverId=NULL)
                            )
              )
       ),
                    
# Panel 2: Kinome Screens
       tabPanel("Kinome Screens",
               
                bsCollapse(multiple = F, id = "kinome", open="kinome1",
                           bsCollapsePanel("Collapse", id='kinome1',
                                fluidRow( 
                                        column(width=3, 
                                             h5('1. Select Samples'),
                                             selectInput('kinome_selected_samples', NULL, choices=unique(kinomeData$condition),
                                                        selected = unique(kinomeData$condition)[1], selectize=T, multiple=T),
                                             br(),
                                             tags$a(href=global_cellLines_metadata_link,target="_blank", "cell line metadata")
                                          ),
                                        column(width=3,
                                               h5('2. Select Kinase Family'),
                                               selectInput('kinome_selected_kinaseFamily', NULL, choices = c('ALL',unique(kinomeData$Family)),
                                                           selected = 'ALL', selectize=T, multiple=T)
                                        ),
                                        column(width=3,
                                               h5('3. Select Genes'),
                                               selectInput('kinome_selected_genes', NULL, choices = c(unique(kinomeData$Gene)),
                                                           selectize=T, multiple=T)
                                        )
                                )
                           )
                ),
                hr(),
                #showOutput("kinome_barPlot", "nvd3")
                #downloadButton(outputId = "download_kinomeBarPlot", label="download plot"),
                
                fluidRow(column(width=6,
                                plotOutput("kinome_barPlot_1")
                                ),
                         column(width=3,
                                h4('Filter by:'),
                                br(),
                                sliderInput("kinome_var_threshold", "a. % variation", 
                                            min=0, max=200, value=c(0,200), step=5),
                                plotOutput("kinome_var_histogram"),
                                
                                hr(),
                                sliderInput('kinome_ratio_excludeRegion', 'c. Exclude Ratio',
                                            min=1, max=2, value=c(1), step=0.1),
                                #plotOutput("kinome_ratio_histogram"),
                                
                                hr(),
                                br(),
                                h5('test: experiment plot'),
                                ggvisOutput('kinome_peptide_scatter_plot_1')
                                
                         ),
                         column(width=3,
                                br(),br(),br(),
                                sliderInput("kinome_uniqPeptideCount_threshold", "b. Unique Peptide Counts", 
                                            min=0, max=50, value=1, step=1),
                                plotOutput("kinome_uniqPeptides_histogram"),
                                hr()
                                
                                
                         )
                )
       ),

      #Tab panel 3
      tabPanel("Public Data",
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
                    tags$textarea(id="custom_gene_list",rows=8,style="width:300px"),
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
     footer = list(img(src="synodos-banner.jpg", height="60", width="120"),
                    helpText('Please report any bugs and/or user experience based feedback at apratap@sagebase.org')
                   )
   
  )#END nav bar
) #END ShinyUI
