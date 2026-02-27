#-----------------
#ATR_FTIR
#-----------------
library(mdatools)
library(pmp)       # tylko PQN
library (Rcpm)   # Inna wersja PQN

setwd ("c:/Users/ar123/Dropbox/Kamila/Dane_uporzadkowane/IR")   # AR dom
#setwd ("c:/Users/ar/Dropbox/Kamila/Dane_uporzadkowane/IR")      # AR praca
#setwd ("c:/Users/ar/Dropbox/Kamila/Macierze/ATR")
#
# Zakres pełny: 10-1750 (566.9692-3922.501) cm-1
#
dolne=as.integer(10)      # 10 (566.9692)   Zakres pełny
gorne=as.integer(1750)      # 1750 (3922.501)
chart_name="Full range"
#chart_name="Raw spectra"
#dolne=as.integer(253)      # 253 (1035.587)   Zakres obciety
#gorne=as.integer(316)      # 316 (1157.08)
dolne=as.integer(235)      # 253 (1002.803)   Zakres obciety Molecules
gorne=as.integer(391)      # 316 (1301.715)
chart_name="Fingerprint range"
size_xy=1.
cex_labels=0.8
line_odstep=2.4
cex_labels=1.

df=read.csv2("macierz_ziola ATR_BLOST.csv", head=TRUE, dec = ".")
df1=read.csv2("macierz_ziola ATR_BLOST.csv", head=FALSE, dec = ".")

cat ("Wybrany zakres: dolne,gorne: ",df1[1,dolne]," ",df1[1,gorne],"\n")

przed_obrobka = df[,dolne:gorne]
labels=as.numeric(df1[1,dolne:gorne])
cl1=factor(df[,1791])
attr(przed_obrobka, "xaxis.values") = labels

po_obrobce = przed_obrobka
przed_obrobka_smoothed = prep.savgol(przed_obrobka, width = 15, porder = 1)
po_obrobce_snv = prep.norm(przed_obrobka,"snv")
po_obrobce_pqn = prep.norm(przed_obrobka,"pqn")
po_obrobce_pqn_v1=pqn_normalisation(df=przed_obrobka,classes=labels,qc_label='all')
X=as.matrix(przed_obrobka)
X = as.matrix (przed_obrobka_smoothed)  # Wygładzanie
po_obrobce_pqn_v2_Rcpm=pqn(X, n = "median", QC = NULL)   # To jest z RCPM
po_obrobce_smoothed_der = prep.savgol(po_obrobce_pqn_v2_Rcpm, width = 15, porder=2, dorder = 1)

par(mfrow = c(2, 2))
#par(mfrow = c(1, 1))
#par(mfrow = c(1, 2))
#
#   Fingerprint
#
yl=substitute(paste(bold('Absorbance')))
ee = mdaplot(przed_obrobka, type = "l", xlab="", main = "Fingerprint range ",cgroup=cl1,lab.col = "black",
             lab.cex=cex_labels,yticks = c(20),yticklabels = c("Y"),
             xticks = c(20),xticklabels = c("A"))
axis(1,cex.axis=size_xy)
axis(2,cex.axis=size_xy)
xl <- (bquote('Wavenumber ['~ cm^-1 ~"]"))   # Bez bolda
xl <- (bquote(bold('Wavenumber ['~ cm^-1 ~"]")))   # Z boldem
mtext(xl,side=1, line=line_odstep, cex=size_xy)
mtext(yl,side=2, line=line_odstep, cex=size_xy)
#
#   Smoothed
#
yl=substitute(paste(bold('Absorbance')))
ee = mdaplot(przed_obrobka_smoothed, type = "l", xlab="", main = "Smoothed",cgroup=cl1,lab.col = "black",
             lab.cex=cex_labels,yticks = c(20),yticklabels = c("Y"),
             xticks = c(20),xticklabels = c("A"))
axis(1,cex.axis=size_xy)
axis(2,cex.axis=size_xy)
xl <- (bquote('Wavenumber ['~ cm^-1 ~"]"))   # Bez bolda
xl <- (bquote(bold('Wavenumber ['~ cm^-1 ~"]")))   # Z boldem
mtext(xl,side=1, line=line_odstep, cex=size_xy)
mtext(yl,side=2, line=line_odstep, cex=size_xy)
#
#   PQN
#
yl=substitute(paste(bold('Absorbance')))
attr(po_obrobce_pqn_v2_Rcpm, "xaxis.values") = labels
mdaplot(po_obrobce_pqn_v2_Rcpm, type = "l", xlab="", main = "PQN",cgroup=cl1,lab.col = "black",
             lab.cex=cex_labels,yticks = c(20),yticklabels = c("Y"),
             xticks = c(200000),xticklabels = c("A"))
axis(1,cex.axis=size_xy)
axis(2,cex.axis=size_xy)
xl <- (bquote('Wavenumber ['~ cm^-1 ~"]"))   # Bez bolda
xl <- (bquote(bold('Wavenumber ['~ cm^-1 ~"]")))   # Z boldem
mtext(xl,side=1, line=line_odstep, cex=size_xy)
mtext(yl,side=2, line=line_odstep, cex=size_xy)
#
#   PQN derivative                           
#
yl=substitute(paste(bold('Absorbance der.')))
attr(po_obrobce_smoothed_der, "xaxis.values") = labels
mdaplot(po_obrobce_pqn_v2_Rcpm, type = "l", xlab="", main = "PQN derivative",cgroup=cl1,lab.col = "black",
        lab.cex=cex_labels,yticks = c(20),yticklabels = c("Y"),
        xticks = c(200000),xticklabels = c("A"))
axis(1,cex.axis=size_xy)
axis(2,cex.axis=size_xy)
xl <- (bquote('Wavenumber ['~ cm^-1 ~"]"))   # Bez bolda
xl <- (bquote(bold('Wavenumber ['~ cm^-1 ~"]")))   # Z boldem
mtext(xl,side=1, line=line_odstep, cex=size_xy)
mtext(yl,side=2, line=line_odstep, cex=size_xy)
