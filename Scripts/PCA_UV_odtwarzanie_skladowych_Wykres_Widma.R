#-----------------
#UV-VSIS
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
par(mfrow = c(2, 2))
#
# Jeden z dwóch - SNV lub PQN
#
po_obrobce_snv = prep.snv(przed_obrobka)
po_obrobce_pqn = prep.norm(przed_obrobka,"pqn")

gg=as.integer(gorne-dolne+1)
do_pca_raw = przed_obrobka[1:93,1:gg]     #  Raw data
do_pca_snv = po_obrobce_snv[1:93,1:gg]    #  SNV
do_pca_pqn = po_obrobce_pqn[1:93,1:gg]    #  PQN
X=do_pca_pqn
po_obrobce_pqn_v2_Rcpm=pqn(X, n = "median", QC = NULL)   # To jest z RCPM
# do_pca_wygl = prep.savgol(do_pca_pqn, width = 15, porder = 1)
do_pca_wygl = prep.savgol(po_obrobce_pqn_v2_Rcpm, width = 15, porder = 1)
do_pca_poch = prep.savgol(po_obrobce_pqn_v2_Rcpm , width = 15, porder = 2, dorder = 1)
################################################################
#
#  Tu decyzja o tym czy:
# 1. Surowe dane (raw)
# 2. Tylko normalizacja (SNV/PQN)
# 3. PQN + Wygladzane/pochodne (m_wygl_poch)
#
# Poczatek obszaru decyzyjnego
do_pierwszego_wykresu = przed_obrobka
#do_pierwszego_wykresu = przed_obrobka_smoothed
#do_pierwszego_wykresu = po_obrobce_snv
do_pierwszego_wykresu = po_obrobce_pqn_v2_Rcpm
#do_pierwszego_wykresu = po_obrobce_smoothed_der
chart_name="Original PQN"
attr(do_pierwszego_wykresu, "xaxis.values") = labels
mdaplot(do_pierwszego_wykresu, type = "l", ylab = "Absorbance", xlab=bquote ('Wavenumber  ['~ cm^-1 ~"]"), 
        main = chart_name,cgroup=cl1,lab.col = "black",lab.cex=1.)
#
#  PCA
#
ncomp <- as.integer(1)
m_jedyne = pca(do_pierwszego_wykresu, ncomp, scale = FALSE,
               info = "PCA model", center= FALSE)
loadingsT=t(m_jedyne$loadings)
scores=m_jedyne$calres$scores
spectra_odtworzone <- scores %*% loadingsT
numf = as.character(ncomp)
#gdzie_pisac=paste(sciezka_do_outputow,"IR_odtw_PCA_ncomp",numf,".csv",sep="")
#write.csv(spectra_odtworzone,"odtworzone_widmo_IR_PQN2.csv")
#ncc=ncol(spectra_odtworzone)

chart_name="Reproduced ncomp=1"
dfo=spectra_odtworzone
attr(dfo, "xaxis.values") = labels
mdaplot(dfo, type = "l", ylab = "Absorbance", xlab=bquote ('Wavenumber  ['~ cm^-1 ~"]"),
        main = chart_name,cgroup=cl1,lab.col = "black",lab.cex=1.)
#
#   Druga reprodukcja
#
ncomp=2
m_jedyne = pca(do_pierwszego_wykresu, ncomp, scale = FALSE, info = "PCA model", center= FALSE)
loadingsT=t(m_jedyne$loadings)
scores=m_jedyne$calres$scores
spectra_odtworzone <- scores %*% loadingsT
chart_name="Reproduced ncomp=2"
dfo=spectra_odtworzone
attr(dfo, "xaxis.values") = labels
mdaplot(dfo, type = "l", ylab = "Absorbance", xlab=bquote ('Wavenumber  ['~ cm^-1 ~"]"),
        main = chart_name,cgroup=cl1,lab.col = "black",lab.cex=1.)
#
#   Trzecia reprodukcja
#
ncomp=3
m_jedyne = pca(do_pierwszego_wykresu, ncomp, scale = FALSE, info = "PCA model", center= FALSE)
loadingsT=t(m_jedyne$loadings)
scores=m_jedyne$calres$scores
spectra_odtworzone <- scores %*% loadingsT
chart_name="Reproduced ncomp=3"
dfo=spectra_odtworzone
attr(dfo, "xaxis.values") = labels
mdaplot(dfo, type = "l", ylab = "Absorbance", xlab=bquote ('Wavenumber  ['~ cm^-1 ~"]"),
        main = chart_name,cgroup=cl1,lab.col = "black",lab.cex=1.)

