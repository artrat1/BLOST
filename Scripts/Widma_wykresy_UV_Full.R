#-----------------
#ATR_FTIR
#-----------------
library(mdatools)
library(pmp)       # tylko PQN
library (Rcpm)   # Inna wersja PQN

#setwd ("C:/Users/user/dropbox/Kamila/Dane_uporzadkowane/UV")   # Laptop K
try (setwd ("C:/Users/ar123/dropbox/Kamila/Dane_uporzadkowane/UV"))    # Dom AR
try (setwd ("C:/Users/ar/dropbox/Kamila/Dane_uporzadkowane/UV"))      # Praca AR
#
# Zakres pełny:  cm-1
#
      
dolne=as.integer(6)      #  500   Zakres pełny
gorne=as.integer(148)     # 212
chart_name="Full range"
chart_name="Raw spectra"
dolne=as.integer(33)      #     Zakres obciety
gorne=as.integer(99)      #  
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
mdaplot(przed_obrobka, type = "l", ylab = "Absorbance", xlab=bquote ('wavelength [nm]'),
        main = chart_name,cgroup=cl1,lab.col = "black",lab.cex=1.)
mdaplot(przed_obrobka_smoothed, type = "l", ylab = "Absorbance", xlab=bquote ('wavelength [nm]'), 
        main = "Smoothed" ,cgroup=cl1,lab.col = "black",lab.cex=1.)
attr(po_obrobce_pqn_v2_Rcpm, "xaxis.values") = labels
mdaplot(po_obrobce_pqn_v2_Rcpm, type = "l", ylab = "Absorbance",
        xlab=bquote ('wavelength [nm]'), main = "PQN" ,cgroup=cl1,lab.col = "black",lab.cex=1.)
attr(po_obrobce_smoothed_der, "xaxis.values") = labels
mdaplot(po_obrobce_smoothed_der, type = "l", ylab = "Absorbance der.", xlab=bquote ('wavelength [nm]'), main = "PQN derivative" 
        ,cgroup=cl1,lab.col = "black",lab.cex=1.)

gg=as.integer(gorne+1-dolne)
stop()
#*********************************************************************************
#
# Teraz będzie podział na zioła, rysowane będzie to co zdefiniowane poniżej
#
#*********************************************************************************
to_aktualnie_rysujemy = przed_obrobka
to_aktualnie_rysujemy = po_obrobce_smoothed_der
par(mfrow = c(2, 3))
#
# Tutaj widma dla poszczegolnych ziół
#
df_lav = to_aktualnie_rysujemy[1:9,1:gg]
attr(df_lav, "xaxis.name") = "Wave number, cm-1"
attr(df_lav, "xaxis.values") = labels
mdaplot(df_lav, type = "l", ylab = "Absorbance",xlab=bquote ('wavelength [nm]'), 
        main = "Lavender",cgroup=cl1[1:9],lab.col = "black",lab.cex=1.)

df_bas = to_aktualnie_rysujemy[43:63,1:gg]
attr(df_bas, "xaxis.name") = "Wave number, cm-1"
attr(df_bas, "xaxis.values") = labels
mdaplot(df_bas, type = "l", ylab = "Absorbance",xlab=bquote ('wavelength [nm]'), 
        main = "Basil",cgroup=cl1[43:63],lab.col = "black",lab.cex=1.)
#
df_oreg = to_aktualnie_rysujemy[64:84,1:gg]
attr(df_oreg, "xaxis.name") = "Wave number, cm-1"
attr(df_oreg, "xaxis.values") = labels
mdaplot(df_oreg, type = "l", ylab = "Absorbance",xlab=bquote ('wavelength [nm]'), main = "Oregano",
        cgroup=cl1[64:84],lab.col = "black",lab.cex=1.)
#
df_sage = to_aktualnie_rysujemy[10:42,1:gg]
attr(df_sage, "xaxis.name") = "Wave number, cm-1"
attr(df_sage, "xaxis.values") = labels
mdaplot(df_sage, type = "l", ylab = "Absorbance",xlab=bquote ('wavelength [nm]'), 
    main = "Sage",cgroup=cl1[10:42],lab.col = "black",lab.cex=1.)
#
df_thyme = to_aktualnie_rysujemy[85:93,1:gg]
attr(df_thyme, "xaxis.name") = "Wave number, cm-1"
attr(df_thyme, "xaxis.values") = labels
mdaplot(df_thyme, type = "l", ylab = "Absorbance",xlab=bquote ('wavelength [nm]'), main = "Thyme",
        cgroup=cl1[85:93],lab.col = "black",lab.cex=1.)
#
#   Teraz wykres widm uśrednionych - tylko 5 linii każda to uśrednione widmo innego zioła
#
#par(mfrow = c(2, 2))
mm <- matrix(, nrow = 5, ncol = gg)
xx=data.frame(mm)
for (i in 1:gg)  {
  xx[1,i]=mean(df_lav[,i])
  xx[2,i]=mean(df_bas[,i])
  xx[3,i]=mean(df_oreg[,i])
  xx[4,i]=mean(df_sage[,i])
  xx[5,i]=mean(df_thyme[,i])
}
cl2 <- factor(c("B","L","O","S","T"))
attr(xx, "xaxis.values") = labels
#
#  Linia poniżej to tylko do Fig. 1 (linia przekopiowana z góry)
#
#par(mfrow = c(1, 2))
#mdaplot(przed_obrobka, type = "l", ylab = "Absorbance", xlab=bquote ('wavelength [nm]'), main = chart_name,cgroup=cl1)
#
#   A to jest faktycznie średnia, tego nie mozna przenieść wyżej
#
mdaplot(xx, type = "l", ylab = "Absorbance", xlab=bquote ('wavelength [nm]'),
        main = "Spectra averaged",cgroup=cl2,lab.col = "black",lab.cex=1.)

