shinyUI(navbarPage("Synodos NF2 Data Explorer",
                   tabPanel("Drug Screens",
                            drugScreenModuleUI("demo") 
                   )
        )
)