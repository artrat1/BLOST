#-----------------
#ATR_FTIR
#-----------------
library(mdatools)
setwd ("c:/Users/ar123/Dropbox/Kamila/Dane_uporzadkowane/IR")   # AR dom
#setwd ("c:/Users/ar/Dropbox/Kamila/Dane_uporzadkowane/IR")      # AR praca
#setwd ("c:/Users/ar/Dropbox/Kamila/Macierze/ATR")
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

attr(przed_obrobka, "xaxis.name") = "Wavelength, nm"
attr(przed_obrobka, "xaxis.values") = labels

po_obrobce = przed_obrobka
po_obrobce_snv = prep.norm(przed_obrobka,"snv")
po_obrobce_pqn = prep.norm(przed_obrobka,"pqn")
par(mfrow = c(2, 2))

gg=as.integer(gorne+1-dolne)

po_obrobce_df=przed_obrobka                # Tu wybieramy sposób obróbki (bez, SNV, PQN)
#po_obrobce_df = po_obrobce_snv[1:93,1:gg]
po_obrobce_df = po_obrobce_pqn[1:93,1:gg]

attr(po_obrobce_df, "xaxis.name") = "Wavelength, nm"
attr(po_obrobce_df, "xaxis.values") = labels
#Colors<-c("#FDA75F","#F1E19D","#E5AC4D","#FDF07A","#FDB98C")
cl1=factor(df[,117])
mdaplot(przed_obrobka, type = "l", ylab = "Absorbance", main = "Original",cgroup=cl1)
mdaplot(po_obrobce_df, type = "l", main = "after SNV", ylab = "Absorbance",cgroup=cl1)

po_obrobce_df1 = po_obrobce_df[1:93,1:gg]
pspectra =  po_obrobce_df[1:93,1:gg]

#pspectra = prep.savgol(po_obrobce_df, width = 15, porder = 1)          # wygładzanie
dpspectra=pspectra
#po_obrobce_df = pspectra  [1:93,1:as.integer(gg)]
#dpspectra = prep.savgol(pspectra , width = 15, porder = 2, dorder = 1)  # pochodna
po_obrobce_df1 = dpspectra  [1:93,1:as.integer(gg)]

attr(po_obrobce_df, "xaxis.name") = "Wavelength, nm"
attr(po_obrobce_df, "xaxis.values") = labels
attr(po_obrobce_df1, "xaxis.name") = "Wavelength, nm"
attr(po_obrobce_df1, "xaxis.values") = labels

mdaplot(po_obrobce_df, type = "l", main = "Smoothed", ylab = "Absorbance")
mdaplot(po_obrobce_df1, type = "l", main = "Pochodna", ylab = "Absorbance")

idx2=seq(2,93,3)  # Tak było na początku, training set  drugie elementy
# 100%(2):  Zakres: 30-93, PQN , ncomp= 3,6,4,2,2
idx1=seq(1,92,3)  # Pierwsze elementy
# 100%(1): Zakres: 30-93, PQN , ncomp= 3,6,4,2,2
idx3=seq(3,93,3)  # trzecie elementy

idx=idx3    #   Tu jest wybor training setu !!

po_obrobce_dfx=po_obrobce                   # tylko SNV
po_obrobce_dfx=po_obrobce_df               # SNV + wygładzania
po_obrobce_dfx=po_obrobce_df1              # SNv + wygładzanie + pochodna

Xc=po_obrobce_dfx[-idx,1:as.integer(gg)]
Xt=po_obrobce_dfx[idx,1:as.integer(gg)]

cc=factor(po_obrobce_dfx[-idx,1])

ct=factor(po_obrobce_dfx[idx,1])

x.l=Xc[1:6,]
x.s=Xc[7:28,]
x.b=Xc[29:42,]
x.o=Xc[43:56,]
x.t=Xc[57:62,]
y.l=Xt[1:3,]
y.s=Xt[4:14,]
y.b=Xt[15:21,]
y.o=Xt[22:28,]
y.t=Xt[29:31,]

nl=as.integer(3)
ns=as.integer(4)
no=as.integer(4)
nb=as.integer(2)
nt=as.integer(2)

m.l = simca(x.l, "lavender", nl)
m.ls = selectCompNum(m.l, nl)
m.s = simca(x.s, "sage", ns)
m.ss = selectCompNum(m.s, ns)
m.b = simca(x.b, "basil", nb)
m.bs = selectCompNum(m.b, nb)
m.o = simca(x.o, "origanum", no)
m.os = selectCompNum(m.o, no)
m.t = simca(x.t, "thyme", nt)
m.ts = selectCompNum(m.t, nt)

mm=simcam(list(m.ls,m.ss,m.bs,m.os,m.ts))
summary(mm)
res = predict(mm, Xt, ct)
layout(1)
plotPredictions(mm$calres) 

par(mfrow = c(2, 3))
#plotModelDistance(mm, 1)
#plotModelDistance(mm, 2)
#plotModelDistance(mm, 3)
#plotModelDistance(mm, 4)
#plotModelDistance(mm, 5)

