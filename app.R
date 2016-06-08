source("global.R")

ui <- fluidPage(
  drugScreenModuleUI("demo") 
)

server <- function(input,output,session){
  #test <- callModule(testModule,"demo")
  callModule(drugScreenModule,id = "demo",session = session, summarizedData = summarizedData, rawData = rawData)#,test)
}

shinyApp(ui,server)