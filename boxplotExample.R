myData<-sizeAndExpressionData[,]#969 x 56514
myData<-sizeAndExpressionData[,-(1:2)]
myData[1:10,1:10]
myData<-t(myData)
mySizes<-sizeAndExpressionData$stage_event_tnm_categories
mySizes<-as.factor(mySizes)
boxplot(myData["ENSG00000062822",] ~mySizes)
myData["ENSG00000062822",]
length(myData["ENSG00000062822",])
length(mySizes)
boxplot(log2(myData["ENSG00000062822",]) ~mySizes)
