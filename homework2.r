#MED 2017Z Lab4, godz. 16
#autor: Michał Sypetkowski, Krzysztof Kukiełka

library(arulesSequences)
#https://archive.ics.uci.edu/ml/datasets/Diabetes
#pobranie danych
download.file('http://staff.ii.pw.edu.pl/~gprotazi/dydaktyka/dane/diab_trans.data','diab_trans.data')
#wczytanie danych do ramki danych
diab.df <- read.csv("diab_trans.data", header=TRUE, stringsAsFactors = FALSE)

#każdą ze zmiennych ciągłych poddano dyskretyzacji dzieląc ją na przedziały
#dla pomiarów które analizujemy ważne jest czy wyniki znajdują się w określonych zakresach 
#przykładowo norma poziomu glukozy we krwi przed śniadaniem jest określona jako <100, natomiast wartości z zakresu 70-130 
#mogą już wskazywać na cukrzycę u danej osoby. dlatego wszystkie wartości poniżej 70 można zagregować i zastąpić wartością dyskretną np. "br1"
#zmienne podzielono na przedziały tak, aby umożliwiały rozróżnienie poszczególnych stanów - 'norma', 'podwyższony' określonych przez Joslin Diabetes Center
#Dzięki temu będzie możliwe skostruowanie reguł sekwencyjnych w oparciu o przedziały, a nie dokładne wartości liczbowe.
mask <- (diab.df$code=='id_33')
insData = diab.df[mask,]
insData$value <- cut(insData$value, c(0,6,7, +Inf), labels = c("lowIns", "medIns", "highIns"))

mask <- (diab.df$code=='id_58')
bbreakfastData = diab.df[mask,]
bbreakfastData$value <- cut(bbreakfastData$value, c(0,70,100,130, +Inf), labels = c("br1", "br2", "br3", "br4"))

mask <- (diab.df$code=='id_62')
supperData = diab.df[mask,]
supperData$value <- cut(supperData$value, c(0,70,110,130, +Inf), labels = c("lun1", "lun2", "lun3", "lun4"))

#z uwagi na to, że zdarzenia o id=48 nie miały określonych progów norm, porównano histogramy podobnych 
#zmiennych i wybrano taką samą dyskretyzację jak dla zdarzeń o id=62
hist(diab.df[diab.df$code=='id_48',]$value)
hist(diab.df[diab.df$code=='id_62',]$value)
hist(diab.df[diab.df$code=='id_58',]$value)
mask <- (diab.df$code=='id_48')
unknownData = diab.df[mask,]
unknownData$value <- cut(unknownData$value, c(0,70,110,130, +Inf), labels = c("unk1", "unk2", "unk3", "unk4"))

#zdarzenia o id=34 miały inny rozkład wartości niż pozostałe zmienne, dlatego zdecydowano o doborze przedziałów na podstawie histogramu zmiennej
#wstępna analiza wykazała występowanie tzw. outlier-ów, czyli wartości odstających od większości. Były to niewielka grupa obserwacji o wartości >50
#Na podstawie histogramu dla wartości <50 dobrano następujące podziały przedziałów: 10, 25, 40 
print(unique(diab.df[which(diab.df$code=='id_34' & diab.df$value),]$value))
hist(diab.df[which(diab.df$code=='id_34' & diab.df$value<50),]$value)
mask <- (diab.df$code=='id_34')
nphData = diab.df[mask,]
nphData$value <- cut(nphData$value, c(0,10,25,40, +Inf), labels = c("nph1", "nph2", "nph3", "nph4"))

#po dyskretyzacji każdej ze zmiennych dokonano złączenia zbiorów 
diab.df <- rbind(insData, bbreakfastData, supperData, unknownData, nphData)

#jako transakcję traktowano powiązane ze sobą id zdarzenia oraz wartość przedziałową zdarzenia
#w tym celu dokonano połączenia id zdarzenia oraz wartość poprzez sklejenie wartości tekstowych
#np. dla kodu='id_33' i value='highIns' otrzymano codeval='id_33_highIns' po sklejeniu
diab.df$codeval <- paste(diab.df$code, diab.df$value, sep='_')
diab.df <- diab.df[, c("patient_id", "time_sek", "codeval")]

View(diab.df)

#ostatecznie tabela bazowa do wykrywania sekwencji zawiera trzy kolumny: patient_id, time_sek oraz codeval
#tabela została posortowana rosnąco po patient_id i time_sek (na potrzeby algorytmu cspade)
diab.df <- diab.df[with(diab.df, order(patient_id, time_sek)),]

#zapis tabeli do pliku 
write.table(diab.df, "diab_trans2.data", sep = "," , row.names = FALSE, col.names = FALSE )

#wczytanie danych w postaci transkacji 
#?read_baskets
diabSeq <- read_baskets(con = "diab_trans2.data", sep =",", info = c("sequenceID","eventID"))
View(as(diabSeq,"data.frame"))

#ustawienie parametrów
# w 1 iteracji wykonujemy algorytm z domyślnymi parametrami
seqParam = new ("SPparameter", support = 0.5, maxsize = 4, mingap=600, maxgap =172800, maxlen = 3 )
patSeq= cspade(diabSeq,seqParam, control = list(verbose = TRUE, tidLists = FALSE, summary= TRUE))

#odkrycie reguł
seqRules = ruleInduction(patSeq,confidence = 0.8)

length(seqRules)
#podsumowanie 
summary(seqRules)
#prezentacja przykĹadowych reguĹ
seqRules <- sort(seqRules, decreasing = TRUE, by="lift")
inspect(head(seqRules,10))

# w powyższej iteracji uzyskano 5 reguł sekwencyjnych o poniższych parametrach:
# <{"id_34_nph1"},                              
# {"id_62_lun2"}>    => <{"id_34_nph1"}>    0.500000  1.0000000 1.523810 
# 2 <{"id_34_nph1"},                              
# {"id_34_nph1"}>    => <{"id_34_nph1"}>    0.593750  1.0000000 1.523810 
# 3 <{"id_34_nph1"},                              
# {"id_58_br2"}>     => <{"id_34_nph1"}>    0.500000  0.9696970 1.477633 
# 4 <{"id_33_highIns"},                           
# {"id_33_highIns"}> => <{"id_33_highIns"}> 0.515625  0.9428571 1.471777 
# 5 <{"id_34_nph1"},                              
# {"id_33_lowIns"}>  => <{"id_34_nph1"}>    0.593750  0.9500000 1.447619 
# są to reguły o najwyższych wartościach liftu z. 
# narzucony w pierwszym kroku parametr maxlen może być zbyt mały - w dalszych krokach szukano dłuższych sekwencji, które mają określone następniki z 
# dużą pewnością. Wtedy można się spodziewać lepszych rezultatów prognostyki. Krótkie sekwencje mogą nieść większe ryzyko przypadkowości
# W kolejnym kroku zmniejszono także wartość parametru maxgap do 18 godzin, tak aby wyszukiwać bliższe sobie zdarzenia - pozwoli to uniknąć wyszukiwania
# sekwencji, które są powtórzeniem tych samych zdarzeń w kolejnych dniach

# 2 iteracja
seqParam = new ("SPparameter", support = 0.5, maxsize = 4, mingap=600, maxgap =64800, maxlen = 5 )
patSeq= cspade(diabSeq,seqParam, control = list(verbose = TRUE, tidLists = FALSE, summary= TRUE))
seqRules = ruleInduction(patSeq,confidence = 0.8)
seqRules <- sort(seqRules, decreasing = TRUE, by="lift")
inspect(head(seqRules,20))

# 1 <{"id_34_nph1"},                              
# {"id_33_lowIns"},                            
# {"id_58_br4"}>     => <{"id_34_nph1"}>    0.531250  1.0000000 1.523810 
# 2 <{"id_34_nph1"},                              
# {"id_58_br4"}>     => <{"id_34_nph1"}>    0.562500  1.0000000 1.523810 
# 3 <{"id_33_lowIns"},                            
# {"id_34_nph1"},                              
# {"id_58_br4"}>     => <{"id_34_nph1"}>    0.546875  1.0000000 1.523810 
# 4 <{"id_34_nph1"},                              
# {"id_34_nph1"},                              
# {"id_58_br4"}>     => <{"id_34_nph1"}>    0.500000  1.0000000 1.523810 
# 5 <{"id_34_nph1"},                              
# {"id_33_lowIns"},                            
# {"id_34_nph1"},                              
# {"id_33_lowIns"}>  => <{"id_34_nph1"}>    0.531250  1.0000000 1.523810
# widać, że w najlepszych wynikach pojawiły się reguły, które mają taki sam następnik. Może wynikać to z przyjętego progu wsparcia.
#Przyjęty próg wsparcia 0.5 oznacza, że reguła pojawia się dla połowy pacjentów. W następnej iteracji zwiększono próg wsparcia, aby
# wyszukiwać reguły, które występują u większej części pacjentów - wtedy z większą szansą będzie można powiedzieć, że reguła, która nie 
# zadziałała u danego pacjenta może być anomalią

# 3 iteracja
seqParam = new ("SPparameter", support = 0.8, maxsize = 4, mingap=600, maxgap =64800, maxlen = 5 )
patSeq= cspade(diabSeq,seqParam, control = list(verbose = TRUE, tidLists = FALSE, summary= TRUE))
seqRules = ruleInduction(patSeq,confidence = 0.8)
seqRules <- sort(seqRules, decreasing = TRUE, by="lift")
inspect(head(seqRules,20))

# 1 <{"id_58_br4"},                             
# {"id_33_lowIns"}> => <{"id_58_br4"}>    0.843750  1.0000000 1.0847458 
# 2 <{"id_33_lowIns"},                          
# {"id_58_br4"}>    => <{"id_62_lun4"}>   0.828125  0.9636364 1.0453005 
# 3 <{"id_58_br4"}>    => <{"id_62_lun4"}>   0.875000  0.9491525 1.0295892 
# 4 <{"id_62_lun4"}>   => <{"id_58_br4"}>    0.875000  0.9491525 1.0295892 
# 5 <{"id_33_lowIns"},                          
# {"id_58_br4"}>    => <{"id_33_lowIns"}> 0.828125  0.9636364 1.0278788 
# w wynikach pojawiły się reguły krótsze niż wcześniej. Być może przyjęty support 0.8 jest zbyt silnym ograniczeniem.
# poza tym pogorszyły się także wyniki liftu z 1.5 do ok 1.05, co oznacza pogorszenie jakości predykcyjnej 
# w następnej iteracji został zmniejszony próg liftu