drugScreenModuleUI <- function(id){
  ns <- NS(id)
  tagList(
    myHeader <- dashboardHeader(title="Drug Screen", disable=TRUE),
    mySidebar <- dashboardSidebar(disable=TRUE),
    myBody <- dashboardBody(
      tags$head(tags$style(HTML('
                                .btn {
                                float:right;
                                }'))),
      fluidRow(
        box(width=12, status='primary', collapsible=TRUE, 
            collapsed=FALSE, solidHeader=TRUE,
            title = tagList(shiny::icon("th-list", lib="glyphicon"),
                            "Data Selection"),
            column(width = 5,
                   h4('1. Select Samples'), 
                   selectInput(ns('samples'),NULL, choices = unique(summarizedData$sample),
                               selectize=T, multiple=T,selected = unique(summarizedData$sample)[1:3])
            ),
            column(width = 7,
                   h4('2. Select Drugs'),
                   h5("Select by drug name"),
                   selectInput(ns('selected_drugs'),NULL, choices = unique(summarizedData$drug),
                               selectize=T, multiple=T, selected = unique(summarizedData$drug)[1:3]),
                   uiOutput(ns("target_class"))
            ),
            actionButton(ns("updateButton"), "Update")
        ),
        box(width = 12, status = 'warning', collapsible = TRUE,
            collapsed = TRUE, solidHeader = TRUE,
            title = tagList(shiny::icon("filter", lib="glyphicon"),
                            "Filter"),
            #uiOutput(ns("filters"))
            column(width=3,
                   # Max Response filter 
                   sliderInput(ns('maxR_filter'), 'Max Response', 
                               min=1, max=100, value=c(1, 100), 
                               step=10, round=TRUE)
            ),
            column(width=3,
                   # IC50 filter
                   sliderInput(ns('ic50_filter'), 'IC50 (uM)',
                               min=10, max=100, value=c(1, 10), 
                               step=5)
            ),
            column(width=3,
                   # AC50 filter
                   uiOutput(ns("filter1"))
            ),
            column(width=3,
                   # Curve class filter
                   uiOutput(ns("filter2"))
            )
        )
      ),
      
      fluidRow(
        tabBox(width = 12,
               tabPanel("Max Response",
                        plotOutput(ns("drug_max_resp"))
               ),
               tabPanel("IC50",
                        plotOutput(ns("drugScreen_IC50_plot"))
               ),
               tabPanel("AC50",
                        plotOutput(ns("drugScreen_AC50_plot"))
               ),
               tabPanel("Dose Response",
                        helpText("If more than 8 drugs are selected, only the first 8 drugs will be showing."),
                        #checkboxInput(ns("replicate"), "use \"replicate\""),
                        #br(),
                        plotOutput(ns("doseResp_plot"))
               ),
               tabPanel("Data",
                        downloadButton(ns("downloadData")),
                        br(),
                        br(),
                        br(),
                        dataTableOutput(ns("drugScreen_dataTable"))
               ),
               tabPanel("QC",
                        h5("Density Histograms Plots"),
                        plotOutput(ns("QC_plots"),height="500px")
               )
        )
      )
      )
    )
  dashboardPage(header=myHeader, sidebar=mySidebar, body=myBody,
                skin = "blue")
} 

shinyUI(navbarPage("Synodos Data Explorer",
             tabPanel("Drug Screens",
                      drugScreenModuleUI("demo") 
             )
        )
)