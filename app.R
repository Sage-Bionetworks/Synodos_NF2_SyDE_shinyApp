source("global.R")

ui <- fluidPage(
  drugScreenModuleUI("demo") 
)

server <- function(input,output,session){
  #test <- callModule(testModule,"demo")
  callModule(drugScreenModule,id = "demo",session = session, drugScreenData)#,test)
}

shinyApp(ui,server)