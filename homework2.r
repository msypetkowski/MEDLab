#MED 2017Z Lab4, godz. 16
#autor: Michał Sypetkowski, Krzysztof Kukiełka

library(arulesSequences)
#https://archive.ics.uci.edu/ml/datasets/Diabetes
#pobranie danych
download.file('http://staff.ii.pw.edu.pl/~gprotazi/dydaktyka/dane/diab_trans.data','diab_trans.data')
#wczytanie danych do ramki danych
diab.df <- read.csv("diab_trans.data", header=TRUE, stringsAsFactors = FALSE)


mask <- (diab.df$code=='id_33')
insData = diab.df[mask,]
insData$value <- cut(insData$value, c(0,6,7, +Inf), labels = c("lowIns", "medIns", "highIns"))

mask <- (diab.df$code=='id_58')
bbreakfastData = diab.df[mask,]
bbreakfastData$value <- cut(bbreakfastData$value, c(0,70,100,130, +Inf), labels = c("br1", "br2", "br3", "br4"))

mask <- (diab.df$code=='id_62')
supperData = diab.df[mask,]
supperData$value <- cut(supperData$value, c(0,70,110,130, +Inf), labels = c("lun1", "lun2", "lun3", "lun4"))

diab.df <- rbind(insData, bbreakfastData, supperData)

diab.df$codeval <- paste(diab.df$code, diab.df$value, sep='_')
diab.df <- diab.df[, c("patient_id", "time_sek", "codeval")]

View(diab.df)

diab.df <- diab.df[with(diab.df, order(patient_id, time_sek)),]

#przykĹad zapisu danych do pliku - usuniÄcie wiersza nagĹĂłwka
write.table(diab.df, "diab_trans2.data", sep = "," , row.names = FALSE, col.names = FALSE )

#wczytanie danych w postaci transkacji 
?read_baskets
diabSeq <- read_baskets(con = "diab_trans2.data", sep =",", info = c("sequenceID","eventID"))
View(as(diabSeq,"data.frame"))

#ustawienie parametrĂłw
seqParam = new ("SPparameter",support = 0.5, maxsize = 4, mingap=600, maxgap =172800, maxlen = 3 )
patSeq= cspade(diabSeq,seqParam, control = list(verbose = TRUE, tidLists = FALSE, summary= TRUE))

#odkrycie reguĹ
seqRules = ruleInduction(patSeq,confidence = 0.8)

length(seqRules)
#podsumowanie 
summary(seqRules)
#prezentacja przykĹadowych reguĹ
inspect(head(seqRules,100))
