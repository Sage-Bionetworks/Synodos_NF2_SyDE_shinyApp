shinyServer(function(input, output, session) {
  session$sendCustomMessage(type="readCookie", message=list())
  foo <- observeEvent(input$cookie, {
    
    synLogin(sessionToken=input$cookie)
    message(sprintf("Welcome, %s", synGetUserProfile()$userName))
    source("getData.R")

    callModule(drugScreenModule,id = "demo",session = session, summarizedData = summarizedData, rawData = rawData,tag="demo")
    callModule(expressionViewerModule,id = "RNA",session = session,RNAseq,pathways_list,tag="RNA")
    callModule(kinomeViewerModule,id = "kinome",session = session,Syn5.Syn1.base=Syn5.Syn1.base,HS01.HS11.base=HS01.HS11.base,kinometx,tag="kinome")
  })
})
  
