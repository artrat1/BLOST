library (mdatools)
#-----------------
#ATR_FTIR
#-----------------
try(setwd ("c:/Users/ar123/Dropbox/Kamila/Dane_uporzadkowane/IR"))   # AR dom
try(setwd ("c:/Users/ar/Dropbox/Kamila/Dane_uporzadkowane/IR") )     # AR praca

#df=read.csv2("BLOST_980-1200.csv")
dolne=as.integer(30)
gorne=as.integer(93)

df=read.csv2("BLOST_980-1200.csv", head=TRUE, dec = ".")
df1=read.csv2("BLOST_980-1200.csv", head=FALSE, dec = ".")
#df=read.csv2("macierz_ziola ATR_BLOST.csv", head=TRUE, dec = ".")
#df1=read.csv2("macierz_ziola ATR_BLOST.csv", head=FALSE, dec = ".")

cat ("Wybrany zakres: dolne,gorne: ",df1[1,dolne]," ",df1[1,gorne])

przed_obrobka = df[,dolne:gorne]
labels=as.numeric(df1[1,dolne:gorne])

po_obrobce = przed_obrobka
po_obrobce_snv = prep.norm(przed_obrobka,"snv")
po_obrobce_pqn = prep.norm(przed_obrobka,"pqn")
par(mfrow = c(2, 2))

gg=as.integer(gorne+1-dolne)

wygladzone_pqn = prep.savgol(po_obrobce_pqn, width = 15, porder = 2)  # wygladzone
pochodna_pqn   = prep.savgol(wygladzone_pqn, width = 15, porder = 2, dorder = 1)  # pochodna


po_obrobce_dfx=przed_obrobka                   # raw
po_obrobce_dfx=po_obrobce_snv               # SNV 
po_obrobce_dfx=po_obrobce_pqn              # PQN
#po_obrobce_dfx=pochodna_pqn                # PQN + SG
#

gg=as.integer(gorne+1-dolne)
k=0.3
kk=1./k
idx = round(seq(1, nrow(po_obrobce_dfx), by = kk))

xc=po_obrobce_dfx[-idx,1:as.integer(gg)]
xt=po_obrobce_dfx[idx,1:as.integer(gg)]
cc=factor(df[-idx,117])
ct=factor(df[idx,117])
cc_all=factor(df[,117])

#train <- xt
#test <- xc
train <- xc
test <- xt
# Załóżmy, że 'train' to ramka danych zawierająca dane treningowe
# 'test' to ramka danych zawierająca dane testowe

# Indeksy kolumn i klas
feature_cols <- 2:116
class_col <- 117

########################################################################################
#########################################################################################
########################################################################################
############################### KONIEC CZĘSCI WSPOLNEJ #####################################################
#######################################################################################
#
# Wczytanie danych IRIS1
iris1 = data.frame(Species = cc_all,po_obrobce_dfx)
iris2=iris1[,-1]
#
#  IRIS1 zdefiniowany !!!! (Iris2 - Iris1 bez pierwszej kolumny (zmienna FACTOR - BLOST))
#
num_columns=ncol(iris1)
ncc=num_columns+1
training_set=iris2[idx,-ncc]
test_set=iris2[-idx,-ncc]
ct_train=factor(iris1[idx,1])
cc_test=factor(iris1[-idx,1])

x_b <- training_set[ct_train=="B",]
x_l <- training_set[ct_train=="L",]
x_o <- training_set[ct_train=="O",]
x_s <- training_set[ct_train=="S",]
x_t <- training_set[ct_train=="T",]


##############################################################################
#    Tu koniec modyfikacji - stad dalsza praca !!!!!!!
###########################################################################

ncb=as.integer(4)
ncl=as.integer(2)
nco=as.integer(4)
ncs=as.integer(7)
nct=as.integer(1)

m.b = simca(x_b, "B", ncb)
#m.b = selectCompNum(m.set, ncb)
m.l = simca(x_l, "L", ncl)
m.o = simca(x_o, "O", nco)
m.s = simca(x_s, "S", ncs)
m.t = simca(x_t, "T", nct)



mm = simcam(list(m.b, m.l, m.o, m.s, m.t))
summary(mm)

res = predict(mm, test_set, cc_test)
# res = predict(mm, xt, ct)
plotPredictions(res)
summary(res)
ave_accuracy = 1 - res$misclassified[4]
ave_accuracy

