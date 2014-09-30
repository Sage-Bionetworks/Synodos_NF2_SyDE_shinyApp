

#x <- getURL("https://raw.github.com/patilv/rChartsTutorials/master/findata.csv",.opts=list(followlocation=T))
#findata <- read.csv(text=x)
#head(findata)
#head(ddply(melt(findata),.(variable),transform,rescale=rescale(value)))
# These are data regarding NCAA athletic department expenses at public universities. Please see the blog post where these charts were originally used 
# regarding more details on the origins of these data.: http://analyticsandvisualization.blogspot.com/2013/10/subsidies-revenues-and-expenses-of-ncaa.html
#findata=findata[,-c(1:2)] # removing first dummy column - the csv quirk - and second column on Rank.
#findatamelt=ddply(melt(findata),.(variable),transform,rescale=rescale(value))



