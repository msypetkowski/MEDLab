#MED 2017Z Lab4, godz. 16
#autor: Michał Sypetkowski


# 1. Przygotowanie danych / preprocessing =========================================

# download.file('http://staff.ii.pw.edu.pl/~gprotazi/dydaktyka/dane/cardioto_noClass_corr.csv','cardioto_noClass_corr.csv')
ctg_noClass <- read.csv("cardioto_noClass_corr.csv",row.names = 1)

# download.file('http://staff.ii.pw.edu.pl/~gprotazi/dydaktyka/dane/cardioto_all_corr.csv','cardioto_all_corr.csv')
ctg <- read.csv("cardioto_all_corr.csv",row.names = 1)

# Usunięcie elementów z klasą 0 (są tylko 3)
sum(ctg$CLASS == 0) # 3
ctg_noClass <- ctg_noClass[ctg$CLASS!=0,]
ctg <- ctg[ctg$CLASS!=0,]



# 2. Eksperymenty =========================================

getPredictions <- function(data, data.kmeans) {
	t <- table(data$CLASS,data.kmeans$cluster)
	t <- as.data.frame.matrix(t) 
	mapping <- apply(t, 2, which.max)
	predictedClass <- mapping[data.kmeans$cluster]
	stopifnot(length(data$CLASS) == length(predictedClass))
	predictedClass
}

getAccuracy <- function(data, data.kmeans) { 
	predictedClass <- getPredictions(data, data.kmeans)
	correct <- data$CLASS == predictedClass
	sum(correct) / length(correct)
}

doExperiment <- function(data, data_noClass) {
    colors <- c("red", "green", "blue")
    for (a in 1:3) {
        acc <- vector(mode = "double" ,length = 15)
        # wss <- vector(mode = "integer" ,length = 15)

        for (i in 1:15) {
            acc[i] <- 0
            repCount <- 50
            for (rep in 1:repCount) {
                kmeans.group <- kmeans(data_noClass, centers = i, iter.max=a*5)
                # wss[i] <- kmeans.group$tot.withinss
                acc[i] <- acc[i] + getAccuracy(data, kmeans.group)
            }
            acc[i] <- acc[i]/repCount
        }

        # CaĹkowita suma odlegĹoĹci wewnÄtrzgrupowych per liczba grup
        # plot(1:15, wss, type = "b",
        #     xlab = "Liczba grup",
        #     ylab = "CaĹkowita suma odlegĹoĹci wewnÄtrzgrupowych")

        if (a == 1) {
            plot(1:15, acc, type = "b",
                xlab = "Liczba grup",
                ylab = "Accuracy", col=colors[a])
        } else {
            lines(1:15, acc, type = "b",
                xlab = "Liczba grup",
                ylab = "Accuracy", col=colors[a])
        }
    }
    legend(8, 0.35, legend=c("iter.max = 10", "iter.max = 20", "iter.max = 30"),
        col=colors, lty=1:2, cex=1.0)
}

doExperiment(ctg, ctg_noClass)




# ctg.kmeans = kmeans(ctg_noClass,2000, iter.max = 20)
# head(getPredictions(ctg, ctg.kmeans))
# head(ctg$CLASS)
# getAccuracy(ctg, ctg.kmeans)
# 
# 
# # podgląd na 2 pierwszych atrybutach
# plot(ctg_noClass[, 1:2], col = ctg$CLASS)
# plot(ctg_noClass[1:2], col = getPredictions(ctg, ctg.kmeans))
# plot(ctg_noClass[1:2], col = getPredictions(ctg, ctg.kmeans) == ctg$CLASS)
