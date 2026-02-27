#-----------------
#Uv-Vis PCA
#-----------------
library(mdatools)
#rstudioapi::writeRStudioPreference("data_viewer_max_columns", 200L)
#setwd ("C:/Users/user/Documents/z dropboxa")

#setwd ("C:/Users/ar/dropbox/Kamila/Dane_uporzadkowane/UV")       # Praca AR
setwd ("C:/Users/ar123/dropbox/Kamila/Dane_uporzadkowane/UV")       # Dom AR
#
#  Poniżej - obcięcia zakresu
#
#  Full range (230-500)
#
dolne=as.integer(4)    #  Sensowne od 4 (ponizej zmienne nienumeryczne)
gorne=as.integer(139)  #  Sensowne do 139 (powyzej szumy)
#
# 370-442
#
#dolne=as.integer(33)   # Sensowne od 4 (ponizej zmienne nienumeryczne)
#gorne=as.integer(69)   # Sensowne do 139 (powyzej szumy)
#
# 390-462
#
dolne=as.integer(23)   # Sensowne od 4 (ponizej zmienne nienumeryczne)
gorne=as.integer(59)   # Sensowne do 139 (powyzej szumy)
#
# 248-500
#dolne=as.integer(4)     # 500
#gorne=as.integer(130)   # 248
#
df=read.csv2("dataUV.csv", head=TRUE, dec = ",")
df1=read.csv2("dataUV.csv", head=FALSE, dec = ",")
przed_obrobka = df[1:93,dolne:gorne]
labels=as.numeric(df1[1,dolne:gorne])

po_obrobce_snv = prep.snv(przed_obrobka)
po_obrobce_pqn = prep.norm(przed_obrobka,"pqn")

mdaplot(przed_obrobka, type = "l", main = "Original")
mdaplot(po_obrobce_snv, type = "l", main = "after SNV")
mdaplot(po_obrobce_pqn, type = "l", main = "after PQN")

gg=as.integer(gorne-dolne+1)
do_pca_raw = przed_obrobka[1:93,1:gg]     #  Raw data
#
do_pca_snv = po_obrobce_snv[1:93,1:gg]    #  SNV
do_pca_pqn = po_obrobce_pqn[1:93,1:gg]    #  PQN
# 

do_pca_wygl = prep.savgol(do_pca_pqn, width = 15, porder = 1)
do_pca_poch = prep.savgol(do_pca_wygl , width = 15, porder = 2, dorder = 1)

m_raw = pca(do_pca_raw, 7, scale = TRUE, info = "PCA model")
m_snv_only=pca(do_pca_snv, 7, scale = TRUE, info = "PCA model")
m_pqn_only = pca(do_pca_pqn, 7, scale = TRUE, info = "PCA model")
m_wygl_poch=pca(do_pca_poch, 7, scale = TRUE, info = "PCA model")
#
#  Tu decyzja o tym czy wygladzane/pochodne (m_wygl_poch) czy nie (m_raw)
#
# Poczatek obszaru decyzyjnego
m=m_raw        # Dane surowe bez pbrobki
m=m_snv_only   # tylko SNV
m=m_pqn_only  # PQN, ale bez wygladzania i pochodnych
m=m_wygl_poch     # z wygladzaniem i pochodnymi
#  Koniec obszaru decyzyjnego
################################################################

print(m)
par(mfrow = c(1, 1))
mdaplot(m$res$cal$scores, type = "p", show.labels = TRUE, show.lines = c(0, 0))

g1 <- factor(df[, "plant"])
#g2 <- factor(people[, "Region"], labels = c("S", "M"))
#
#  Plot
#
par(mfrow = c(2, 2))
p= plotScores(m$res$cal, c(1,2),show.labels = FALSE, cgroup=g1)
plotConfidenceEllipse(p)

p= plotScores(m$res$cal, c(1,3),show.labels = FALSE, cgroup=g1)
plotConfidenceEllipse(p)

p= plotScores(m$res$cal, c(1,4),show.labels = FALSE, cgroup=g1)
plotConfidenceEllipse(p)

p= plotScores(m$res$cal, c(2,3),show.labels = FALSE, cgroup=g1)
plotConfidenceEllipse(p)

