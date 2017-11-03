shinyServer(function(input, output, session) {
  source("global.R")
  callModule(drugScreenModule,id = "demo",session = session, summarizedData = summarizedData, rawData = rawData,tag="demo")
  callModule(expressionViewerModule,id = "RNA",session = session,RNAseq,pathways_list,tag="RNA")
  callModule(kinomeViewerModule,id = "kinome",session = session,Syn5.Syn1.base=Syn5.Syn1.base,HS01.HS11.base=HS01.HS11.base,kinometx,tag="kinome")
  })