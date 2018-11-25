#MED 2017Z Lab4, godz. 16
#autor: Michał Sypetkowski, Krzysztof Kukiełka

library(arulesSequences)
#https://archive.ics.uci.edu/ml/datasets/Diabetes
#pobranie danych
download.file('http://staff.ii.pw.edu.pl/~gprotazi/dydaktyka/dane/diab_trans.data','diab_trans.data')
#wczytanie danych do ramki danych
diab.df <- read.csv("diab_trans.data", header=TRUE, stringsAsFactors = FALSE)

# 1. Najlepsze reguły =====================================================

# Szukamy reguł sekwencyjnych, które mogą być stosowane w praktyce przez chorych na cukrzycę.
# Dobra reguła musi mieć wysoki współczynnik podniesienia (np. > 2),
# oraz wysokie zafuanie (np > 90%).

# Reguły muszą być krótkie, oraz czytelne / możliwe do zrozumienia przez chorego
# (chociaż intuicyjnie można przypuszczać że długie reguły o wysokim wsp. podniesienia,
# zaufaniu oraz wsparciu (np. > 0.25) nie istnieją w tym zbiorze danych)

# 2. Pre-processing =====================================================

# usuwamy pomiary z wartościami NA
diab.df <- diab.df[!is.na(diab.df$value),]

# Wybieramy te pomiary co występują dostatecznie często
table(diab.df$code)
# pomijamy id występujące mniej niż 1000 razy, oraz te bez wyjaśnienia czym są 
# (https://archive.ics.uci.edu/ml/datasets/diabetes)

# pozostają nam pomiary
# 33 = Regular insulin dose
# 34 = NPH insulin dose
# 35 = UltraLente insulin dose
# 48 = Unspecified blood glucose measurement
# 58 = Pre-breakfast blood glucose measurement
# 60 = Pre-lunch blood glucose measurement
# 62 = Pre-supper blood glucose measurement


#każdą ze zmiennych ciągłych poddano dyskretyzacji dzieląc ją na przedziały
#dla pomiarów które analizujemy ważne jest czy wyniki znajdują się w określonych zakresach 
#przykładowo norma poziomu glukozy we krwi przed śniadaniem jest określona jako <100, natomiast wartości z zakresu 70-130 
#mogą już wskazywać na cukrzycę u danej osoby. dlatego wszystkie wartości poniżej 70 można zagregować i zastąpić wartością dyskretną np. "br1"
#zmienne podzielono na przedziały tak, aby umożliwiały rozróżnienie poszczególnych stanów - 'norma', 'podwyższony' określonych przez Joslin Diabetes Center
#Dzięki temu będzie możliwe skostruowanie reguł sekwencyjnych w oparciu o przedziały, a nie dokładne wartości liczbowe.
mask <- (diab.df$code=='id_33')
data33 = diab.df[mask,]
data33$value <- cut(data33$value, c(-Inf,6,7, +Inf), labels = c("low", "med", "high"))

mask <- (diab.df$code=='id_58')
data58 = diab.df[mask,]
data58$value <- cut(data58$value, c(-Inf,70,100,130, +Inf), labels = c("low", "mlow", "mhigh", "high"))

mask <- (diab.df$code=='id_60')
data60 = diab.df[mask,]
data60$value <- cut(data60$value, c(-Inf,70,110,130, +Inf), labels = c("low", "mlow", "mhigh", "high"))

mask <- (diab.df$code=='id_62')
data62 = diab.df[mask,]
data62$value <- cut(data62$value, c(-Inf,70,110,130, +Inf), labels = c("low", "mlow", "mhigh", "high"))

#z uwagi na to, że zdarzenia o id=48 nie miały określonych progów norm 
# (https://www.joslin.org/info/goals_for_blood_glucose_control.html), porównano histogramy podobnych 
#zmiennych i wybrano taką samą dyskretyzację jak dla zdarzeń o id=62
hist(diab.df[diab.df$code=='id_48',]$value)
hist(diab.df[diab.df$code=='id_62',]$value)
hist(diab.df[diab.df$code=='id_58',]$value)
mask <- (diab.df$code=='id_48')
data48 = diab.df[mask,]
data48$value <- cut(data48$value, c(-Inf,70,110,130, +Inf), labels = c("low", "mlow", "mhigh", "high"))


#zdarzenia o id=34 miały inny rozkład wartości niż pozostałe zmienne, dlatego zdecydowano o doborze przedziałów na podstawie histogramu zmiennej
#wstępna analiza wykazała występowanie tzw. outlier-ów, czyli wartości odstających od większości. Były to niewielka grupa obserwacji o wartości >50
#Na podstawie histogramu dla wartości <50 dobrano następujące podziały przedziałów: 10, 25, 40 
print(unique(diab.df[which(diab.df$code=='id_34' & diab.df$value),]$value))
hist(diab.df[which(diab.df$code=='id_34' & diab.df$value<50),]$value)
mask <- (diab.df$code=='id_34')
data34 = diab.df[mask,]
data34$value <- cut(data34$value, c(-Inf,10,25,40, +Inf), labels = c("low", "mlow", "mhigh", "high"))

#po dyskretyzacji każdej ze zmiennych dokonano złączenia zbiorów 
diab.df <- rbind(data33, data34, data48, data58, data60, data62)

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


# 3. Eksperymenty =======================================================


#ustawienie parametrów
# szukamy krótkich reguł ze wsparciem conajmniej 50%, dopuszczamy wiele itemów w jednym elemencie sekwencji
# przyjmujemy pomiary w odstępach 20 minutowych jako itemy w tym samym elemencie,
# następników szukamy w zakresie 1 doby.
seqParam = new ("SPparameter", support = 0.5, maxsize=20, mingap=20 * 60, maxgap = 24 * 60 * 60, maxlen = 3 )
?cspade
?ruleInduction
patSeq= cspade(diabSeq,seqParam, control = list(verbose = TRUE, tidLists = FALSE, summary= TRUE))

seqRules = ruleInduction(patSeq,confidence = 0.8)
length(seqRules)
summary(seqRules)
seqRules <- sort(seqRules, decreasing = TRUE, by="lift")
inspect(head(seqRules,10))

# najlepsze (w sensie wsp. podniesienia) 6 reguł to:
# 1 <{"id_34_low"},                         
# {"id_34_low"}>  => <{"id_34_low"}>  0.578125  0.9736842 1.483709 
# 2 <{"id_34_low"},                         
# {"id_58_high"}> => <{"id_34_low"}>  0.562500  0.9729730 1.482625 
# 3 <{"id_33_low"},                         
# {"id_34_low"}>  => <{"id_34_low"}>  0.593750  0.9500000 1.447619 
# 4 <{"id_34_low"},                         
# {"id_33_low"}>  => <{"id_34_low"}>  0.593750  0.9500000 1.447619 
# 5 <{"id_58_high"},                        
# {"id_34_low"}>  => <{"id_34_low"}>  0.546875  0.9459459 1.441441 
# 6 <{"id_58_mlow"},                        
# {"id_34_low"}>  => <{"id_34_low"}>  0.500000  0.9411765 1.434174 

# powyżej uzyskane reguły nie są interesujące. Informują o tym że pomiary id_34 i id_33 są powtarzalne.
# Oznacza to, że chorzy często przyjmują stałe dawki insuliny - co nie jest zaskakujące.
# Towarzyszące okazyjnie występujące itemy jak id_58 to najprawdopodobniej pomiary poprzedzające dawkę insuliny
# (więc chorzy postępują zgodnie z regułą - podwyższony cukier jest wskazaniem do przyjęcia dawki insuliny)
# W regule 6 poziom cukru jest mały, więc przyjmowane dawki mogą być prawdopodobnie zapobiegawcze (np. przed posiłkiem).

# uwzględnijmy tylko te reguły, co nie mają w rhs informacji o przyjmowanych dawkach insuliny
patSeqFiltered = subset(seqRules, ! rhs(seqRules) %in% c('"id_33_low"','"id_33_med"','"id_33_high"', '"id_34_low"'))
patSeqFiltered <- sort(patSeqFiltered, decreasing = TRUE, by="lift")
inspect(head(patSeqFiltered,10))

# 1 <{"id_33_low"},                         
# {"id_58_high"}> => <{"id_62_high"}> 0.859375  0.9821429 1.065375 
# 2 <{"id_58_high"},                        
# {"id_33_low"}>  => <{"id_58_high"}> 0.843750  0.9818182 1.065023 
# 3 <{"id_62_high"},                        
# {"id_33_low"}>  => <{"id_58_high"}> 0.812500  0.9811321 1.064279 
# 4 <{"id_33_low"},                         
# {"id_60_high"}> => <{"id_58_high"}> 0.781250  0.9803922 1.063476 
# 5 <{"id_60_high"},                        
# {"id_33_low"}>  => <{"id_58_high"}> 0.781250  0.9803922 1.063476 
# 6 <{"id_62_high"},                        
# {"id_60_high"}> => <{"id_58_high"}> 0.718750  0.9787234 1.061666 

# otrzymane reguły mają niski współczynnik podniesienia

# szukamy reguł kóre implikują wysoki poziom glukozy we kri
patSeqFiltered = subset(seqRules, rhs(seqRules) %in% c('"id_58_high"','"id_60_high"','"id_62_high"'))
patSeqFiltered <- sort(patSeqFiltered, decreasing = TRUE, by="lift")
inspect(head(patSeqFiltered,10))

# otrzymane reguły również mają (podobnie) niski współczynnik podniesienia



# dopuśćmy dłuższe sekwencje
seqParam = new ("SPparameter", support = 0.5, maxsize = 20, mingap=20 * 60, maxgap = 24 * 60 * 60, maxlen = 5 )
patSeq= cspade(diabSeq,seqParam, control = list(verbose = TRUE, tidLists = FALSE, summary= TRUE))
seqRules = ruleInduction(patSeq,confidence = 0.8)
seqRules <- sort(seqRules, decreasing = TRUE, by="lift")
inspect(head(seqRules,20))

# lhs                 rhs               support confidence    lift 
# 1 <{"id_34_low"},                         
# {"id_34_low"},                         
# {"id_58_high"}> => <{"id_34_low"}>  0.562500          1 1.52381 
# 2 <{"id_58_high"},                        
# {"id_34_low"},                         
# {"id_34_low"},                         
# {"id_58_high"}> => <{"id_34_low"}>  0.546875          1 1.52381 
# 3 <{"id_58_high"},                        
# {"id_34_low"},                         
# {"id_33_low"},                         
# {"id_34_low"}>  => <{"id_34_low"}>  0.546875          1 1.52381 
# 4 <{"id_34_low"},                         
# {"id_62_high"},                        
# {"id_33_low"},                         
# {"id_34_low"}>  => <{"id_34_low"}>  0.500000          1 1.52381 
# 5 <{"id_33_low"},                         
# {"id_33_low"},                         
# {"id_34_low"},                         
# {"id_34_low"}>  => <{"id_34_low"}>  0.562500          1 1.52381 

# Podobnie jak poprzednio, zobaczmy reguły bez informacji o dawce insuliny w następniku

patSeqFiltered = subset(seqRules, ! rhs(seqRules) %in% c('"id_33_low"','"id_33_med"','"id_33_high"', '"id_34_low"'))
patSeqFiltered <- sort(patSeqFiltered, decreasing = TRUE, by="lift")
inspect(head(patSeqFiltered,10))

# Otrzymujemy 1 regułę o wsp podniesienia > 1.1
# lhs                 rhs               support confidence     lift 
# 1 <{"id_60_mlow"},                        
# {"id_62_high"},                        
# {"id_33_low"}>  => <{"id_60_high"}> 0.546875  0.9459459 1.121121 

#  ponieważ elementy sekwencji mogą występować w odległości czasowej do 1 doby, reguła jest mało informatywna
# (poza tym ma dość mały współczynnik podniesienia).


# zmniejszmy maxgap do 10h
seqParam = new ("SPparameter", support = 0.5, maxsize = 20, mingap=20 * 60, maxgap = 10 * 60 * 60, maxlen = 5 )
patSeq= cspade(diabSeq,seqParam, control = list(verbose = TRUE, tidLists = FALSE, summary= TRUE))
seqRules = ruleInduction(patSeq,confidence = 0.8)
seqRules <- sort(seqRules, decreasing = TRUE, by="lift")
inspect(head(seqRules,20))


# lhs                 rhs               support confidence      lift 
# 1 <{"id_58_high"}> => <{"id_60_high"}> 0.796875  0.8644068 1.0244821 
# 2 <{"id_33_low"}>  => <{"id_33_low"}>  0.875000  0.9333333 0.9955556 
# 3 <{"id_34_low"}>  => <{"id_33_low"}>  0.609375  0.9285714 0.9904762 
# 4 <{"id_60_low"}>  => <{"id_62_high"}> 0.718750  0.9019608 0.9783981 
# 5 <{"id_62_high"},                        
# {"id_33_low"}>  => <{"id_58_high"}> 0.515625  0.8918919 0.9674760 

# 58 = Pre-breakfast blood glucose measurement
# 60 = Pre-lunch blood glucose measurement
# 62 = Pre-supper blood glucose measurement

# 1 reguła ma lift > 1:
# Gdy ktoś miał wysoki poziom cukru przed śniadaniem,
# to prawdopodobieństwo że będzie miał również wysoki poziom cukru przed obiadem jest zwiększone


# jak do tąd nie uzyskaliśmy ciekawych reguł dla minimalnego wsparcia 0.5

# dopuśćmy minimalne wsparcie 0.25
seqParam = new ("SPparameter", support = 0.25, maxsize = 20, mingap=20 * 60, maxgap = 10 * 60 * 60, maxlen = 5 )
patSeq= cspade(diabSeq,seqParam, control = list(verbose = TRUE, tidLists = FALSE, summary= TRUE))
seqRules = ruleInduction(patSeq,confidence = 0.8)
seqRules <- sort(seqRules, decreasing = TRUE, by="lift")
inspect(head(seqRules,20))

# lhs                 rhs               support confidence     lift 
# 1 <{"id_48_high"},                        
# {"id_58_high"},                        
# {"id_60_high"},                        
# {"id_62_high"}> => <{"id_48_high"}> 0.250000  1.0000000 3.368421 
# 2 <{"id_48_high"},                        
# {"id_58_high"},                        
# {"id_60_high"}> => <{"id_48_high"}> 0.250000  1.0000000 3.368421 
# 3 <{"id_62_high"},                        
# {"id_48_high"},                        
# {"id_58_high"},                        
# {"id_60_high"}> => <{"id_48_high"}> 0.250000  1.0000000 3.368421 
# 4 <{"id_48_high"},                        
# {"id_58_high"},                        
# {"id_62_high"}> => <{"id_48_high"}> 0.265625  0.9444444 3.181287 
# 5 <{"id_34_low"},                         
# {"id_33_low"},                         
# {"id_33_low"},                         
# {"id_62_mlow"}> => <{"id_34_low"}>  0.265625  1.0000000 1.523810 
# 6 <{"id_34_low"},                         
# {"id_58_high"},                        
# {"id_33_low"},                         
# {"id_62_high"}> => <{"id_34_low"}>  0.265625  1.0000000 1.523810 

# pierwsze 5 reguł ma wysoki wsp. podniesienia

# 48 = Unspecified blood glucose measurement
# 58 = Pre-breakfast blood glucose measurement
# 60 = Pre-lunch blood glucose measurement
# 62 = Pre-supper blood glucose measurement

# W poprzednikach reguł, mamy sekwencje oznaczające wysoki poziom cukru podczas doby
# (id 58, 60, 62 - przed śniadaniem, obiadem i kolacją).
# Wysoki cukier utrzymujący się przez dłuższy czas (około doby), stanowi dobrą przesłankę (jest przyczyną)
# że następny pomiar (wieczorem, prawdopodobnie po kolacji lub po obiedzie -
# ostatni element sekwencji lhs to najczęściej 60 lub 62) również będzie wysoki.


# dopuśćmy bardzo długie długie reguły
seqParam = new ("SPparameter", support = 0.25, maxsize = 20, mingap=20 * 60, maxgap = 10 * 60 * 60, maxlen = 10 )
patSeq= cspade(diabSeq,seqParam, control = list(verbose = TRUE, tidLists = FALSE, summary= TRUE))
seqRules = ruleInduction(patSeq,confidence = 0.8)
seqRules <- sort(seqRules, decreasing = TRUE, by="lift")
inspect(head(seqRules,20))

# nie uzyskujemy reguł z lift > 2 z długością > 6 (intuicja z punktu 1. się sprawdza).


# 4. Wnioski ==============================================

# Wnioski z odpowiednich ekperymentów są zawarte w sekcji 3.

# Ogólne obserwacje:

# Prawdopodobnie nie istnieją użyteczne reguły o wysokim wsparciu ( > 0.5), z wysokim współczynnikiem podniesienia.
