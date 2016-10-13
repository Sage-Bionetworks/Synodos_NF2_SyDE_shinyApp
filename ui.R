shinyUI(navbarPage("Synodos NF2 Data Explorer",
	#google analytics
	header=list(tags$head(includeScript("www/google_analytics.js"))),
   	tabPanel("Drug Screens",
            drugScreenModuleUI(id = "demo", data = summarizedData) 
    	), 
	tabPanel("NCATS Single Agent Drug Screen",
	         drugScreenModuleUI(id = "ncats", data = ncatsData) 
	),
	tabPanel("NCATS Combination Drug Screen",
		combinationDrugScreenModuleUI(id = "ncats_comb",combinedData))	
  )#END navnarPage
)#END shinyUI

