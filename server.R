shinyServer(function(input, output, session) {
  callModule(drugScreenModule,id = "demo",session = session, summarizedData = summarizedData, rawData = rawData,tag="demo")
  callModule(expressionViewerModule,id = "RNA",session = session,RNAseq,pathways_list,tag="RNA")
  callModule(kinomeViewerModule,id = "kinome",session = session,kinome,tag="kinome")
  })