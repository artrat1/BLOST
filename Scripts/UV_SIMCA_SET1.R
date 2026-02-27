library ("mdatools")
# setwd ("C:/Users/user/dropbox/Kamila/Dane_uporzadkowane/UV")   # Laptop K
try (setwd ("C:/Users/ar123/dropbox/Kamila/Dane_uporzadkowane/UV"))  # Dom AR
try (setwd ("C:/Users/ar/dropbox/Kamila/Dane_uporzadkowane/UV"))      # Praca AR
size_xy=1.

line_odstep=2.4
cex_labels=0.9
dolne=as.integer(33)
gorne=as.integer(99)

df=read.csv2("dataUV.csv", head=TRUE, dec = ",")
df1=read.csv2("dataUV.csv", head=FALSE, dec = ",")
przed_obrobka = df[,dolne:gorne]
labels=as.numeric(df1[1,dolne:gorne])

po_obrobce_snv = prep.norm(przed_obrobka,"snv")
po_obrobce_pqn = prep.norm(przed_obrobka,"pqn")

gg=as.integer(gorne+1-dolne)

przed_wygladzeniem =  przed_obrobka[1:93,1:gg]        #  bez obrobki
przed_wygladzeniem =  po_obrobce_snv[1:93,1:gg]    #  SNV
przed_wygladzeniem =  po_obrobce_pqn[1:93,1:gg]    #  PQN

wygladzone = prep.savgol(przed_wygladzeniem, width = 15, porder = 2)          # wygładzanie
pochodna = prep.savgol(wygladzone , width = 15, porder = 2, dorder = 1)  # pochodna

idx2=seq(2,93,3)  # Tak było na początku, training set  drugie elementy
# 100%(2): Zakres: 33-99, bez obróbki, ncomp=5,4,3,2,1
# OPT(2): Zakres: 33-99, PQN, ncomp=6,4,6,3,1
idx1=seq(1,92,3)  # Pierwsze elementy
# 100%(1): Zakres: 33-93, bez obróbki,       ncomp(SOBLT)=8,3,5,2,1    
# 100%(2): Zakres: 33-93, PQN oraz PQNSG,    ncomp(SOBLT)=10,4,8,3,1
# (0.968 - Sage, reszta 100%) PQN,           ncomp(SOBLT)=5,3,5,3,1
idx3=seq(3,93,3)  # trzecie elementy
# OPT: (100%):Zakres(33-93), bez obrobki, ncomp(SOBLT): 8,3,5,2,1
# OPT: (Sage 0.875, reszta 100%):Zakres(33-93), bez obrobki, ncomp(SOBLT): 4,3,5,2,1
# OPT(PQN) (Sage 0.968, reszta 100%):Zakres(33-99), PQN, ncomp(SOBLT): 5,3,5,3,1
####################################################
idx=idx1    #   Tu jest wybor training setu !!
###################################################
#
#po_obrobce_dfx=przed_obrobka              # Raw
#po_obrobce_dfx = po_obrobce_snv          # tylko SNV
po_obrobce_dfx = po_obrobce_pqn          # tylko PQN
#po_obrobce_dfx = wygladzone              # SNV(PQN) + wygładzania
#po_obrobce_dfx = pochodna                # SNv(PQN) + wygładzanie + pochodna

Xc=po_obrobce_dfx[-idx,1:as.integer(gg)]
Xt=po_obrobce_dfx[idx,1:as.integer(gg)]

cc=factor(df[-idx,2])

ct=factor(df[idx,2])

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

ns=as.integer(6)
no=as.integer(3)
nb=as.integer(5)
nl=as.integer(3)
nt=as.integer(1)

m.b = simca(x.b, "B", nb)
m.bs = selectCompNum(m.b, nb)
m.l = simca(x.l, "L", nl)
m.ls = selectCompNum(m.l, nl)

m.o = simca(x.o, "O", no)
m.os = selectCompNum(m.o, no)
m.s = simca(x.s, "S", ns)
m.ss = selectCompNum(m.s, ns)

m.t = simca(x.t, "T", nt)
m.ts = selectCompNum(m.t, nt)

mm=simcam(list(m.bs,m.ls,m.os,m.ss,m.ts))
summary(mm)
res = predict(mm, Xt, ct)
layout(1)
#plotPredictions(mm$calres, main = 'Training set 3')
plotPredictions(mm$calres, main = "",xticks = c(200),xticklabels = c("A"),
                lab.cex=cex_labels,lab.col='black',cex=0.75) 
xl <- (bquote('Objects'))   # Bez bolda
#xl <- (bquote(bold('Wavenumber ['~ cm^-1 ~"]")))   # Z boldem
axis(1,cex.axis=size_xy)
#axis(2,cex.axis=size_xy)
#mtext(xl,side=1, line=line_odstep, cex=cex_labels)
#mtext(yl,side=2, line=line_odstep, cex=cex_labels)
stop()

par(mfrow = c(2, 3))
plotModelDistance(mm, 1)
plotModelDistance(mm, 2)
plotModelDistance(mm, 3)
plotModelDistance(mm, 4)
plotModelDistance(mm, 5)

par(mfrow = c(2, 2))
library ("ggplot2")
library(grid)
plotCooman(mm, c(1, 2), show.legend = FALSE, main = "Lavender - Sage" )
plotCooman(mm, c(1, 3), show.legend = FALSE, main = "Lavender - Basil")
plotCooman(mm, c(1, 4), show.legend = FALSE, main = "Lavender - Oregano")
plotCooman(mm, c(1, 5), show.legend = FALSE, main = "Lavender - Thyme")
ggsave ("Coomans_UV_Lavender.png")
par(mfrow = c(2, 2))
plotCooman(mm, c(2, 3), show.legend = FALSE, main = "Sage - Basil")
plotCooman(mm, c(2, 4), show.legend = FALSE, main = "Sage - Oregano")
plotCooman(mm, c(2, 5), show.legend = FALSE, main = "Sage - Thyme")
plotCooman(mm, c(3, 4), show.legend = FALSE, main = "Basil - Oregano")

par(mfrow = c(1, 2))
plotCooman(mm, c(3, 5), show.legend = FALSE, main = "Basil - Thyme")
plotCooman(mm, c(4, 5), show.legend = FALSE, main = "Oregano - Thyme")
