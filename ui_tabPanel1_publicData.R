#UI.R for tabPanel 1
#Synodos Data Explorer: Tab Panel 1 : exploring public data

#Tabpanel 1
tabPanel("Public Data",
         progressInit(),
         h3("Explore NF2 Public Datasets"),
         sidebarLayout(
           sidebarPanel(
             h4('1. Select Grouped Pathways '),
             selectInput('grouped_pathways', NULL, choices = names(grouped_pathways),
                         selected=c('NF2_related_pathways'), selectize=T, multiple=T),
             
             br(),br(),br(),
             
             h4('2. Select individual pathway/s'),
             selectInput('pathways', NULL, choices = rownames(pathway_enrichment_scores), multiple=TRUE, 
                         width='400px',selectize=T),
             
             br(), br(), br(),
             h4('3. Select based on genes'), 
             helpText("Accepts HUGO gene names. Gene names may be separated by comma, space,line, comma "),
             tags$textarea(id="custom_gene_list",rows=8,cols=400),
             
             br(), br(), br(), br(), br(), br(),
             get_synodos_banner(width="300")
           ),
           
           mainPanel(
             tabsetPanel(id="pubData_pathway_panel", type="tabs",
                         tabPanel("Pathway Enrichment",
                                  plotOutput("pubData_pathway_heatmap",height="700px",width="auto",hoverId=NULL)
                         ),
                         tabPanel("Expression",
                                  h5("expression data based heatmaps..coming soon")                      
                         )
             )
           ),
           fluid=T
         )
)#END TabPanel("Public Data")