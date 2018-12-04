library(precrec)
library(caret)
library(party)
library(rpart)
library(rpart.plot)
library(e1071)
wine = read.table("wine_red.csv", header = TRUE, sep=";", na.strings= "*")

# 1. Cel badań ############################################################

# Celem badań jest wsparcie klienta wybierającego wina w sklepie. Klasyfikator 
# ma za zadanie wskazywać klientowi klasę wina w podziale wina 'bardzo dobrej jakości' 
# oraz 'posotałe'. W ten sposób klient może na podstawie danych z etykiety oszacować czy 
# dane wino z półki bądź sklepu internetowego jest dobrej jakości, 
# a zarazem czy jest warte wysokiej ceny. Klasyfikator może także pomóc w wyszukiwaniu
# 'okazji cenowych', czyli win wysokiej jakości w niskich cenach.

# 2.Preprocessing ############################################################

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



# 3.Eksperymenty ############################################################

formula <- Q2 ~ fixed.acidity + volatile.acidity + citric.acid + residual.sugar + chlorides + free.sulfur.dioxide + total.sulfur.dioxide + density + pH + sulphates + alcohol

getErrorRate <- function(classifier, testData) {
	testPred <- predict(wineCtree, newdata = testData)
	t <- table(testData$Q2,testPred)
	print(t)
	# confusionMatrix(testPred,testData$Q2)
	ret <- (t[2] + t[3]) / (t[1] + t[2] + t[3] + t[4])
	print('Error rate:')
	print(ret)
	ret
}

drawPrecRec <- function(classifier, testData) {
	testPrecProb <- predict(wineCtree, newdata = testData, type="prob")
	testLabels <- cbind(testData$Q2 == 'A', testData$Q2 == 'B')
	testPredDf <- data.frame(matrix(unlist(testPrecProb), nrow=nrow(testLabels), byrow=T))
	colnames(testPredDf) <- c('A', 'B')
	sscurves <- evalmod(scores = testPredDf$A, labels = testData$Q2 == 'A')
	plot(sscurves)
}


# drzewo decyzyjne (ctree) z domyślnymi parametrami
wineCtree <- ctree(formula, data=trainData)
plot(wineCtree)
getErrorRate(wineCtree, testData)
drawPrecRec(wineCtree, testData)
# dostajemy error rate 0.1106472, jak widać na wykresie precission recall,
# dla wysokiego współczynnika odzysku otrzymujemy mniejszą precyzję
# jakość a jest determinowana na podstawie decyzji alcohol > 11.5 i sulphates > 0.68

# spróbujmy ograniczyć minimalną ilość elementów przypodziale
wineCtree <- ctree(formula, data=trainData, controls=ctree_control(minsplit=200))
plot(wineCtree)
getErrorRate(wineCtree, testData)
drawPrecRec(wineCtree, testData)
# otrzymujemy error rate 0.1356994, czyli wynik jest gorszy.
# Ze względu na małą ilość win klasy 'A', potrzebny jest bardziej szczegółowy podział
# (który został uniemożliwiony poprzez ustawienie minsplit=200)

# spróbujmy zmniejszyć ograniczenie minimalnej wartości testu statystycznego potrzebnej do wykonania podziału
wineCtree <- ctree(formula, data=trainData, controls=ctree_control(mincriterion=0.65))
plot(wineCtree)
getErrorRate(wineCtree, testData)
drawPrecRec(wineCtree, testData)
# nie otrzymujemy lepszego error rate
wineCtree <- ctree(formula, data=trainData, controls=ctree_control(mincriterion=0.6))
plot(wineCtree)
getErrorRate(wineCtree, testData)
drawPrecRec(wineCtree, testData)
# error rate zmalał do 0.09394572, pojawił się nowy liść o klasie 'A'

# spróbujmy zmienić kryterium podziału (inny niż "Bonferroni")
wineCtree <- ctree(formula, data=trainData, controls=ctree_control(testtype="Univariate"))
plot(wineCtree)
getErrorRate(wineCtree, testData)
drawPrecRec(wineCtree, testData)
wineCtree <- ctree(formula, data=trainData, controls=ctree_control(testtype="Univariate", mincriterion=0.94))
plot(wineCtree)
getErrorRate(wineCtree, testData)
drawPrecRec(wineCtree, testData)
# nie ma widocznej poprawy

# najlepszy uzyskany klasyfikator jest z rodzajem testu "Bonferroni" i mincriterion=0.6
# (względem domyślnych wartości parametrów)

# 4. Wnioski #############################

# Aby uzyskać precyzyjne pomiary precyzji, należy użyć coross-validation do pomiarów,
# a także przeprowadzić eksperymenty "w pętli" dla różnych wartości różnych parametrów
# (np. iloczyn kartezjański możliwych wartości parametrów, lub algorytm ewolucyjny)
