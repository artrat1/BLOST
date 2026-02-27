#-----------------
#ATR_FTIR
#-----------------
library(mdatools)
#setwd ("c:/Users/ar123/Dropbox/Kamila/Dane_uporzadkowane/IR")   # AR dom
setwd ("c:/Users/ar/Dropbox/Kamila/Dane_uporzadkowane/IR")      # AR praca
#setwd ("c:/Users/ar/Dropbox/Kamila/Macierze/ATR")
#df=read.csv2("BLOST_980-1200.csv")
dolne=as.integer(30)
gorne=as.integer(93)
size_xy=1.
cex_labels=0.8
line_odstep=2.4
cex_labels=1.
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

idx2=seq(2,93,3)  # Tak było na początku, training set  drugie elementy
# 100%(2):  Zakres: 30-93, PQN , ncomp= 3,6,4,2,2
idx1=seq(1,92,3)  # Pierwsze elementy
# 100%(1): Zakres: 30-93, PQN , ncomp= 3,6,4,2,2
idx3=seq(3,93,3)  # trzecie elementy

idx=idx1    #   Tu jest wybor training setu !!

po_obrobce_dfx=przed_obrobka                   # raw
po_obrobce_dfx=po_obrobce_snv               # SNV 
po_obrobce_dfx=po_obrobce_pqn              # PQN
#po_obrobce_dfx=pochodna_pqn                # PQN + SG

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

m.l = simca(x.l, "L", nl)
m.ls = selectCompNum(m.l, nl)
m.s = simca(x.s, "S", ns)
m.ss = selectCompNum(m.s, ns)
m.b = simca(x.b, "B", nb)
m.bs = selectCompNum(m.b, nb)
m.o = simca(x.o, "O", no)
m.os = selectCompNum(m.o, no)
m.t = simca(x.t, "T", nt)
m.ts = selectCompNum(m.t, nt)

mm=simcam(list(m.bs,m.ls,m.os,m.ss,m.ts))
summary(mm)
res = predict(mm, Xt, ct)
layout(1)
#plotPredictions(mm$calres,main = 'Training set 1') 
plotPredictions(mm$calres, main = "",xticks = c(200),xticklabels = c("A"),
                lab.cex=cex_labels,lab.col='black') 
xl <- (bquote('Objects'))   # Bez bolda
#xl <- (bquote(bold('objects)))   # Z boldem
axis(1,cex.axis=size_xy)
stop
