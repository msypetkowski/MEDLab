#MED 2017Z Lab4, godz. 16
#autor: Michał Sypetkowski
library(fossil)
library(gsubfn)



# 1. Przygotowanie danych / preprocessing =========================================

# download.file('http://staff.ii.pw.edu.pl/~gprotazi/dydaktyka/dane/cardioto_noClass_corr.csv','cardioto_noClass_corr.csv')
ctg_noClass <- read.csv("cardioto_noClass_corr.csv",row.names = 1)

# download.file('http://staff.ii.pw.edu.pl/~gprotazi/dydaktyka/dane/cardioto_all_corr.csv','cardioto_all_corr.csv')
ctg <- read.csv("cardioto_all_corr.csv",row.names = 1)

# Usunięcie elementów z klasą 0 (są tylko 3)
sum(ctg$CLASS == 0) # 3
ctg_noClass <- ctg_noClass[ctg$CLASS!=0,]
ctg <- ctg[ctg$CLASS!=0,]

ctgScaled <- sweep(ctg_noClass, 2, apply(ctg_noClass, 2, function(x) {sqrt(sum(x*x)/(length(x) - 1))}) , `/`)
ctgCentered <- sweep(ctg_noClass, 2, apply(ctg_noClass, 2 ,mean), `-`)
ctgScaledCentered <- sweep(ctgCentered, 2, apply(ctgCentered, 2, function(x) {sqrt(sum(x*x)/(length(x) - 1))}), `/`)



# 2. Eksperymenty =========================================

getPredictions <- function(classes, data.kmeans) {
	t <- table(classes,data.kmeans$cluster)
	t <- as.data.frame.matrix(t) 
	mapping <- apply(t, 2, which.max)
	predictedClass <- mapping[data.kmeans$cluster]
	stopifnot(length(classes) == length(predictedClass))
	predictedClass
}

getAccuracy <- function(classes, data.kmeans) { 
	predictedClass <- getPredictions(classes, data.kmeans)
	correct <- classes == predictedClass
	sum(correct) / length(correct)
}

maxGroups <- 15

drawPlot <- function(data, ylab, color, first=F, ylim=NA) {
    if (first==T) {
        if (!is.na(ylim)) {
            plot(1:maxGroups, data, type = "b",
                xlab = "Groups count",
                ylab = ylab, col=color, ylim=ylim)
        } else {
            plot(1:maxGroups, data, type = "b",
                xlab = "Groups count",
                ylab = ylab, col=color)
        }
    } else {
        lines(1:maxGroups, data, type = "b",
            xlab = "Groups count",
            ylab = ylab, col=color)
    }
}

doExperiment <- function(classes, data_noClass, repCount=30, randiRepCount=4) {
    acc <- vector(mode = "double" ,length = maxGroups)
    randi <- vector(mode = "double" ,length = maxGroups)
    wss <- vector(mode = "double" ,length = maxGroups)

    for (i in 1:maxGroups) {
        acc[i] <- 0
        randi[i] <- 0
        wss[i] <- 0
        for (rep in 1:repCount) {
            kmeans.group <- kmeans(data_noClass, centers = i, iter.max=20)
            acc[i] <- acc[i] + getAccuracy(classes, kmeans.group)
            if (rep <= randiRepCount) {
                # Indeks randa ma mniejsze odchylenie (pozatym wolno się liczy), więc robimy mniej powtórzeń.
                randi[i] <- randi[i] + rand.index(classes, kmeans.group$cluster)
            }
            wss[i] <- wss[i] + kmeans.group$tot.withinss
        }
        acc[i] <- acc[i]/repCount
        randi[i] <- randi[i]/randiRepCount
        wss[i] <- wss[i]/repCount
    }

    list(acc, randi, wss)
}

colors <- c("red", "green", "blue", "black", "orange")

# Mierzymy dla różnej ilości grup:

# a) Dokładność w kontekście klasyfikacji -
# Dla każdej uzyskanej grupy wybieramy klasę według kryterium:
# najwięcej elementów w tej grupie należy do tej klasy.
# 
# b) Indeks randa.
#
# c) Sumę odległości par elementów wewnątrz grup

list[acc1, randi1, wss1] <- doExperiment(ctg$CLASS, ctgScaled)
list[acc2, randi2, wss2] <- doExperiment(ctg$CLASS, ctgScaledCentered)
list[acc3, randi3, wss3] <- doExperiment(ctg$CLASS, ctg_noClass)

drawPlot(acc1, color=colors[1], ylab="Accuracy", first=T)
drawPlot(acc2, color=colors[2], ylab="Accuracy")
drawPlot(acc3, color=colors[3], ylab="Accuracy")
legend(8, 0.35, legend=c("Scaled", "Centered and Scaled", "Raw"),
    col=colors, lty=1:2, cex=1.0)

drawPlot(randi1,color=colors[1],  ylab="Rand index", first=T)
drawPlot(randi2,color=colors[2],  ylab="Rand index")
drawPlot(randi3,color=colors[3],  ylab="Rand index")
legend(8, 0.35, legend=c("Scaled", "Centered and Scaled", "Raw"),
    col=colors, lty=1:2, cex=1.0)


# Wartości dla "surowych" danych są bardzo duże (bo wiele kolumn ma duże wartości, stąd duże odległości nawet wewnątrz grup)
# drawPlot(wss3,color=colors[3],  ylab="Internal distances sum", first=T, ylim=c(0, maxVal))
maxVal <- max(max(wss1), max(wss2))
drawPlot(wss1,color=colors[1],  ylab="Internal distances sum", first=T, ylim=c(0, maxVal))
drawPlot(wss2,color=colors[2],  ylab="Internal distances sum")
legend(8, 0.35, legend=c("Scaled", "Centered and Scaled", "Raw"),
    col=colors, lty=1:2, cex=1.0)


# 3. Wnioski

# Najleprze grupowanie uzyskujemy dla ilości grup=15 (maksymalna dopuszczalna z założenia),
# przy skalowaniu danych bez wyśrodkowania (dzieleniu przez mean square root po kolumnach).
