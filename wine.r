#MED 2017Z Lab6, godz. 16
#autor: Michał Sypetkowski, Krzysztof Kukiełka

# kod jest częściowo wzięty z projektu (swojego oczywiście) z innego przedmiotu (MOW)


library(rpart)
library(rpart.plot)
# library(rattle)


# download.file('http://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv', 'wine_red.csv');
wine = read.table("wine_red.csv", header = TRUE, sep=";", na.strings= "*")

# 1. Cel badań ############################################################

# Celem badań jest wsparcie klienta wybierającego wina czerwonego w sklepie. Klasyfikator 
# ma za zadanie wskazywać klientowi klasę wina w podziale wina 'bardzo dobrej jakości' 
# oraz 'posotałe'. W ten sposób klient może na podstawie danych z etykiety oszacować czy 
# dane wino z półki bądź sklepu internetowego jest dobrej jakości, 
# a zarazem czy jest warte wysokiej ceny. Klasyfikator może także pomóc w wyszukiwaniu
# 'okazji cenowych', czyli win wysokiej jakości w niskich cenach.

# 2.Preprocessing ############################################################

# przygotowanie atrybutu do klasyfikacji binarnej (grupowanie)
wine$Q2 = lapply(wine[,12], function (x)
{
  if(x >6)  { 1} # pierwsza klasa (dobra jakość)
  else { 2}  # druga klasa (słaba jakość)
})
wine$Q2 = unlist(wine$Q2);
wine$Q2 = as.factor(wine$Q2)

wineFreq = table(wine$Q2)
wineHist <- barplot(wineFreq, xlab = "Jakosc wina", ylab = "Liczba wystapien wina", ylim = c(0,2400)) 
text(x = wineHist, y = wineFreq, label = wineFreq, pos = 3, cex = 1, col = "blue")




# 3.Eksperymenty ############################################################

# Przeprowadzamy powtarzaną walidację krzyżową drzew decyzyjnych tworzonych z różnymi parametrami 



joinDataframes <- function(listOfDf) {
    df <- do.call("rbind", listOfDf)
    df
}


# współczynnik błędu
getClassifierError <- function(predFun, testset, attrName) {
    predicted <- predFun(testset)
    stopifnot(length(predicted) == nrow(testset))
    sum(predicted!= testset[[attrName]]) / nrow(testset)
}


crossValidation <- function(dataset, trainFunction, errorFunction, partitionsCount) {
    n <- partitionsCount
    nr <- nrow(dataset)
    f <- rep(1:ceiling(n), each=floor(nr/n), length.out=nr)
    parts <- split(dataset, f)
    stopifnot(nrow(dataset) - nrow(joinDataframes(parts)) >= 0)
    stopifnot(nrow(dataset) - nrow(joinDataframes(parts)) < nr/n)

    errSamples <- c()
    for (i in 1:n) {
        train <- joinDataframes(parts[-i])
        test <- parts[[i]]
        stopifnot(nrow(train) + nrow(test) == nrow(dataset))
        fit <- trainFunction(train)
        e <- errorFunction(fit, test)
        errSamples <- c(errSamples, c(e))
    }
    errSamples
}


doExperiment <- function(dataset, trainFun, repTimes) {
    errSamples <- vector()
    set.seed(123)
    for (i in 1:repTimes) {
        dataset <- dataset[sample(nrow(dataset)),]
        errFun <- function(predFun, testset) {
            getClassifierError(predFun, testset, "Q2")
        }
        errSamples <- c(errSamples, crossValidation(dataset, trainFun, errFun, 6))
    }
    errSamples
}


# returns predict function with 1 param - testset
trainDecTreeClassifier <- function(dataset, draw=FALSE, cp=0.01, minbucket=15, maxdepth=15, split="gini") {
    param <- (Q2 ~ fixed.acidity + volatile.acidity + citric.acid + residual.sugar + chlorides
                + free.sulfur.dioxide + total.sulfur.dioxide + density + pH + sulphates + alcohol)
    tr <- rpart(param,data = dataset,method="class",control = rpart.control(
                                    cp=cp,
                                    minsplit = minbucket*2,
                                    minbucket= minbucket,
                                    maxdepth=maxdepth
                             ), parms=list(split=split))
    if (draw) {
        # nie działa na r2d3, chociaż lepiej wygląda
        # fancyRpartPlot(tr)

        plot(tr)
        text(tr)
    }
    predFun <- function(testset) {
        p <- predict(tr, testset)
        ret <- as.integer(p[, 1] < 0.5) + 1
    }
    predFun
}



# train many decision trees
decTreeTest <- function(dataset1, repTimes) {
    # parametry dobierane, aby najlepszy klasyfikator miał możliwie wartości "po środku" przedziałów.
    paramCp = list(0.0001, 0.01, 0.1, 0.3)
    paramMinBucket = list(2, 5, 50)
    paramMaxDepth = list(3, 5, 10)
    paramSplit = list("gini", "information")
    records <- apply(expand.grid(paramCp, paramMinBucket, paramMaxDepth, paramSplit), 1, FUN = function(x) {
        samples <- doExperiment(dataset1, function(ds) {
                    trainDecTreeClassifier(ds, cp=x[[1]], minbucket=x[[2]], maxdepth=x[[3]], split=x[[4]])
        }, repTimes=repTimes)
        data.frame(
                    c(x[[1]]),
                    c(x[[2]]),
                    c(x[[3]]),
                    c(x[[4]]),
                    c(mean(samples)),
                    c(sd(samples)),
                    c(min(samples)),
                    c(max(samples))
        )
    })
    records <- joinDataframes(records)
    records <- setNames(records, c("cp", "minbucket", "maxdepth", "split", "err_mean", "err_sd", "err_min", "err_max"))
    records <- records[with(records, order(err_mean)), ]
    print(records)
    records
}

# more repetitions give more precise measures, but requires more computation time
records <- decTreeTest(wine, repTimes=5) # wykonanie z repTimes=5 trwa około pół minuty
trainDecTreeClassifier(wine, draw=T,
    cp=records$cp[1],
    minbucket=records$minbucket[1],
    maxdepth=records$maxdepth[1],
    split=records$split[1]
)


# 4. Wnioski #############################

# Można zaobserwować, że wina klasy "1" mają z reguły większą zawartość alkoholu (według najlepszego uzyskanego klasyfikatora).
# Między innymi najgorszy uzyskany klasyfikator,
# to drzewo z dużą ilością małych liści (jest zbyt duże żeby je narysować) - następuje overfitting.

# Wyciągnięcie ciekawszych wniosków wniosków wymaga wiedzy dziedzinowej (eksperymenty przeprowadzone zostały (pół)automatycznie).
