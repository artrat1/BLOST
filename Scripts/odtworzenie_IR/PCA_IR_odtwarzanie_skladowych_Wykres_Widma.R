#-----------------
#ATR_FTIR
#-----------------
library(mdatools)
library(pmp)       # tylko PQN
library (Rcpm)   # Inna wersja PQN

try (setwd ("c:/Users/ar123/Dropbox/Kamila/Dane_uporzadkowane/IR"))   # AR dom
try (setwd ("c:/Users/ar/Dropbox/Kamila/Dane_uporzadkowane/IR"))      # AR praca
#
# Zakres pełny: 10-1750 (566.9692-3922.501) cm-1
#
ncomp <- as.integer(1)
dolne=as.integer(10)      # 10 (566.9692)   Zakres pełny
gorne=as.integer(1750)      # 1750 (3922.501)
#chart_name="Full range"
chart_name="Raw spectra ATR-FTIR"
dolne=as.integer(253)      # 253 (1035.587)   Zakres obciety
gorne=as.integer(316)      # 316 (1157.08)
#
#*****************************************************************************
dolne=as.integer(236)      # 236 (1002.803)   Zakres obciety Molecules
gorne=as.integer(391)      # 391 (1301.715)
chart_name="Fingerprint range"
#*****************************************************************************
#
df=read.csv2("macierz_ziola ATR_BLOST.csv", head=TRUE, dec = ".")
df1=read.csv2("macierz_ziola ATR_BLOST.csv", head=FALSE, dec = ".")

cat ("Wybrany zakres: dolne,gorne: ",df1[1,dolne]," ",df1[1,gorne],"\n")

przed_obrobka = df[,dolne:gorne]

labels=as.numeric(df1[1,dolne:gorne])
cl1=factor(df[,1791])

#attr(przed_obrobka, "xaxis.name") = expression("Wave number, cm^{-1})
attr(przed_obrobka, "xaxis.values") = labels


po_obrobce = przed_obrobka
przed_obrobka_smoothed = prep.savgol(przed_obrobka, width = 15, porder = 1)
po_obrobce_snv = prep.norm(przed_obrobka,"snv")
po_obrobce_pqn = prep.norm(przed_obrobka,"pqn")
po_obrobce_pqn_v1=pqn_normalisation(df=przed_obrobka_smoothed,classes=labels,qc_label='all')
#X=as.matrix(przed_obrobka)
X = as.matrix (przed_obrobka_smoothed)  # Wygładzanie
po_obrobce_pqn_v2_Rcpm=pqn(X, n = "median", QC = NULL)   # To jest z RCPM
attr(po_obrobce_pqn_v2_Rcpm, "xaxis.values") = labels
po_obrobce_smoothed_der = prep.savgol(po_obrobce_pqn_v2_Rcpm, width = 15, porder=2, dorder = 1)

par(mfrow = c(2, 2))
#par(mfrow = c(2, 2))
#par(mfrow = c(1, 1))
################################################################
#
#  Tu decyzja o tym czy:
# 1. Surowe dane (raw)
# 2. Tylko normalizacja (SNV/PQN)
# 3. QPN + Wygladzane/pochodne (m_wygl_poch)
#
# Poczatek obszaru decyzyjnego
do_pierwszego_wykresu = przed_obrobka
#do_pierwszego_wykresu = przed_obrobka_smoothed
#do_pierwszego_wykresu = po_obrobce_snv
do_pierwszego_wykresu = po_obrobce_pqn_v2_Rcpm
#do_pierwszego_wykresu = po_obrobce_smoothed_der
chart_name="Original PQN"
mdaplot(do_pierwszego_wykresu, type = "l", ylab = "Absorbance", xlab=bquote ('Wavenumber  ['~ cm^-1 ~"]"), 
        main = chart_name,cgroup=cl1,lab.col = "black",lab.cex=1.)
write.csv(do_pierwszego_wykresu,"do_wykresu_1.csv")

#stop()
#
#  PCA
#
#########################
#  Tu decyzja (do_pca)
###########################
#
m_raw = pca(przed_obrobka, ncomp, scale = TRUE, info = "PCA model", center= FALSE)
m_snv_only=pca(po_obrobce_snv, ncomp, scale = TRUE, info = "PCA model, center= FALSE")
#m_pqn_only = pca(do_pca_pqn, ncomp, scale = TRUE, info = "PCA model", center= FALSE)
write.csv("")
m_jedyne = pca(do_pierwszego_wykresu, ncomp, scale = FALSE, info = "PCA model", center= FALSE)
#m_wygl_poch=pca(po_obrobce_smoothed_der, ncomp, scale = TRUE, info = "PCA model", center= FALSE)

################################################################
#
#  Tu decyzja o tym czy:
# 1. Surowe dane (raw)
# 2. Tylko normalizacja (SNV/PQN)
# 3. QPN + Wygladzane/pochodne (m_wygl_poch)
#
# Poczatek obszaru decyzyjnego
m_pierwszy_wykres = m_raw           # Dane surowe bez pbrobki
m_pierwszy_wykres = m_snv_only      # tylko SNV
m_pierwszy_wykres = m_jedyne      # PQN, ale bez wygladzania i pochodnych
#m=m_wygl_poch     # z wygladzaniem i pochodnymi
#  Koniec obszaru decyzyjnego
################################################################
print(m_jedyne)

loadingsT=t(m_jedyne$loadings)
scores=m_jedyne$calres$scores
spectra_odtworzone <- scores %*% loadingsT

numf = as.character(ncomp)
#gdzie_pisac=paste(sciezka_do_outputow,"IR_odtw_PCA_ncomp",numf,".csv",sep="")
write.csv(spectra_odtworzone,"odtworzone_widmo_IR_PQN2.csv")
ncc=ncol(spectra_odtworzone)

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
