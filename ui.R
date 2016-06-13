shinyUI(navbarPage("Synodos NF2 Data Explorer",
	#google analytics
	header=list(tags$head(includeScript("www/google_analytics.js"))),
   tabPanel("Drug Screens",
            drugScreenModuleUI(id = "demo", data = summarizedData) 
    ) # END tabPanel
  )#END navnarPage
)#END shinyUI

