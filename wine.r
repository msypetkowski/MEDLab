wine = read.table("wine_red.csv", header = TRUE, sep=";", na.strings= "*")


# przygotowanie atrybutu do klasyfikacji binarnej (grupowanie)
wine$Q2 = lapply(wine[,12], function (x)
{
  if(x >6)  { "A"}
  else { "B"}   
})
wine$Q2 = unlist(wine$Q2);
wine$Q2 = as.factor(wine$Q2)

wineFreq = table(wine$Q2)
?barplot
wineHist <- barplot(wineFreq, xlab = "Jakosc wina", ylab = "Liczba wystapien wina", ylim = c(0,2400)) 
text(x = wineHist, y = wineFreq, label = wineFreq, pos = 3, cex = 1, col = "blue")


# podział na zbiór trenujący i testowy
sam <- sample(2, nrow(wine), replace=TRUE, prob=c(0.7, 0.3))
trainData <- wine[sam==1,]
testData <- wine[sam==2,]
idTrainData <- unlist(createDataPartition(wine$Q2,p=0.7))
trainData <-wine[idTrainData,]
testData <-wine[-idTrainData,]
table(trainData$Q2)
table(testData$Q2)



formula <- Q2 ~ fixed.acidity + volatile.acidity + citric.acid + residual.sugar + chlorides + free.sulfur.dioxide + total.sulfur.dioxide + density + pH + sulphates + alcohol
wineCtree <- ctree(formula, data=trainData)
plot(wineCtree)
plot(wineCtree, type="simple")
