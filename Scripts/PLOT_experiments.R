#-----------------
#ATR_FTIR
#-----------------
library(mdatools)
library(pmp)       # tylko PQN
library (Rcpm)   # Inna wersja PQN
library("ggplot2")

try(setwd("C:/Users/user/dropbox/Kamila/Dane_uporzadkowane/UV"))   # Laptop K
try (setwd ("C:/Users/ar123/dropbox/Kamila/Dane_uporzadkowane/UV"))    # Dom AR
try (setwd ("C:/Users/ar/dropbox/Kamila/Dane_uporzadkowane/UV"))      # Praca AR
#
# Zakres pełny:  cm-1
#
      
dolne=as.integer(6)      #  500   Zakres pełny
gorne=as.integer(148)     # 212
chart_name="Full range"
chart_name="Raw spectra"
#dolne=as.integer(33)      #     Zakres obciety
#gorne=as.integer(99)      #  
#chart_name="Fingerprint range"

df=read.csv2("dataUV.csv", head=TRUE, dec = ",")
df1=read.csv2("dataUV.csv", head=FALSE, dec = ",")

cat ("Wybrany zakres: dolne,gorne: ",df1[1,dolne]," ",df1[1,gorne],"\n")

przed_obrobka = df[,dolne:gorne]
labels=as.numeric(df1[1,dolne:gorne])
cl1=factor(df[,2])
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
par(mfrow = c(1, 1))
#par(mfrow = c(1, 2))
yl=substitute(paste(bold('Absorbance')))
ww=bquote ('wavelength [nm]')
xl=substitute(paste(bold('Wavelength [nm]')))
ee = mdaplot(przed_obrobka, type = "l", xlab="", main = "Raw spectra UV-Vis",cgroup=cl1,lab.col = "black",
             lab.cex=1.1,yticks = c(20),yticklabels = c("Y"),
             xticks = c(20),xticklabels = c("A"))
axis(1,cex.axis=1.1)
#mtext(Xl,side=1, line=2.4, cex=1.2)
axis(2,cex.axis=1.1)
mtext(xl,side=1, line=2.7, cex=1.1)
mtext(yl,side=2, line=2.3, cex=1.1)

#ee + mdaplot(opacity = 20,grid.col = "blue")

#getPlotColors(ee, col, opacity, cgroup, colmap)
stop
#
# Ponizej oryginal przed przerobkami
# mdaplot(przed_obrobka, type = "l", ylab = "Absorbance", xlab=bquote ('wavelength [nm]'), main = chart_name,cgroup=cl1)
#
