shinyUI(navbarPage("Synodos NF2 Data Explorer",
	#google analytics
	header=list(tags$head(includeScript("www/iframe_resize.js"))),
   tabPanel("Drug Screens",
            drugScreenModuleUI("demo") 
    ) # END tabPanel
  )#END navnarPage
)#END shinyUI

