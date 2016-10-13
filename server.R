shinyServer(function(input, output, session) {
  callModule(drugScreenModule,id = "demo",session = session, summarizedData = summarizedData, rawData = rawData,tag="demo")
  callModule(drugScreenModule,id = "ncats",session = session, summarizedData = ncatsData,tag="ncats")
  callModule(combinationDrugScreenModule,id = "ncats_comb",session = session,combinedData,tag = "ncats_comb")
})
