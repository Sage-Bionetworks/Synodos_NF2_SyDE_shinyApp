shinyUI(navbarPage("Synodos NF2 Data Explorer",
                   header=list(tags$head(includeScript("www/readCookie.js"))),
   tabPanel("Drug Screens",
             drugScreenModuleUI(id = "demo") 
   ),
	 tabPanel("Gene Expression",
	          expressionViewerModuleUI(id = "RNA")
	          ),
	tabPanel("Kinome Perturbation",
	         kinomeViewerModuleUI(id = "kinome")
	)# END tabPanel
  )#END navnarPage
)#END shinyUI

