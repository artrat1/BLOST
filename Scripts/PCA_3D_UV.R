#-----------------
#UV-Vis
#-----------------
library(mdatools)
setwd ("C:/Users/ar/dropbox/Kamila/Dane_uporzadkowane/UV")       # Praca AR
#setwd ("C:/Users/ar123/dropbox/Kamila/Dane_uporzadkowane/UV")       # Dom AR
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
dolne=as.integer(33)   # Sensowne od 4 (ponizej zmienne nienumeryczne)
gorne=as.integer(69)   # Sensowne do 139 (powyzej szumy)
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

attr(przed_obrobka, "xaxis.name") = "Wavelength, nm"
attr(przed_obrobka, "xaxis.values") = labels


po_obrobce_snv = prep.snv(przed_obrobka)
po_obrobce_pqn = prep.norm(przed_obrobka,"pqn")
gg=as.integer(gorne-dolne+1)
#
do_pca_raw = przed_obrobka[1:93,1:gg]     #  Raw data
do_pca_snv = po_obrobce_snv[1:93,1:gg]    #  SNV
do_pca_pqn = po_obrobce_pqn[1:93,1:gg]    #  PQN
do_pca_wygl = prep.savgol(do_pca_pqn, width = 15, porder = 1)
do_pca_poch = prep.savgol(do_pca_wygl , width = 15, porder = 2, dorder = 1)

#
#  Tu decyzja o tym czy wygladzane/pochodne (m_wygl_poch) czy nie (m_raw)
#
# Poczatek obszaru decyzyjnego
xspectra=do_pca_raw      # Dane surowe bez pbrobki
xspectra=do_pca_snv      # tylko SNV
xspectra=do_pca_pqn      # PQN, ale bez wygladzania i pochodnych
xspectra=do_pca_poch     # z wygladzaniem i pochodnymi
#  Koniec obszaru decyzyjnego
################################################################


library(rgl)
#plot3d(pc$x[,1:3], xlab="Component 1", ylab="Component 2", zlab="Component 3", type="n", box=F, axes=T)
#plot3d(m$res$cal$score[,1:3], xlab="Component 1", ylab="Component 2", zlab="Component 3", type="n", box=F, axes=T)
#decorate3d(xlab = "x", ylab = "y", zlab = "z", 
#           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
#           top = TRUE, aspect = FALSE, expand = 1.03)
#spheres3d(m$res$cal$score[,1:3], radius=0.2, col=rep(c("red","green","black")),cgroup=g1)
#grid3d(c("x", "y+", "z"))
#text3d(m$res$cal$score[,1:3][,1:3], text=rownames(m$res$cal$score[,1:3]), adj=1.3)
library ("pca3d")
dc <- defaultPalettePCA3D(n = NULL, transparent = NULL, d3 = FALSE)
pca <- prcomp( xspectra, scale.= TRUE )
summ <- summary(pca)
g1 <- factor(df[, "plant"])#
#  Bez elips i centroidow
#
pca3d( pca, group= g1 ,show.group.labels = FALSE,axes = "black",radius=1.5,show.plane=FALSE,palette=rep(c("steelblue4","aquamarine4","olivedrab3", "orange", "red" )))
pca2d( pca, group= g1 ,show.group.labels = FALSE,axes = "black",radius=1.5,show.plane=FALSE,palette=rep(c("steelblue4","aquamarine4","olivedrab3", "orange", "red" )))
rgl.snapshot('UV_3dplot.png', fmt = 'png')
#
# Bez centroidow z elipsami
#
pca3d( pca, group= g1 ,show.group.labels = FALSE,axes = "black",radius=1.5,show.plane=FALSE,show.ellipses=TRUE,ellipse.ci=0.75,palette=rep(c("steelblue4","aquamarine4","olivedrab3", "orange", "red" )))
pca2d( pca, group= g1 ,show.group.labels = FALSE,axes = "black",radius=1.5,show.plane=FALSE,show.ellipses=TRUE,ellipse.ci=0.75,palette=rep(c("steelblue4","aquamarine4","olivedrab3", "orange", "red" )))
rgl.snapshot('UV_3D_elipsy.png', fmt = 'png')
#
# Z centroidami i elipsami
#
pca3d( pca, group= g1 ,show.group.labels = FALSE,axes = "black",show.centroids=TRUE,radius=1.5,show.plane=FALSE,show.ellipses=TRUE,ellipse.ci=0.75,palette=rep(c("steelblue4","aquamarine4","olivedrab3", "orange", "red" )))
rgl.snapshot('3dplot_centroidy.png', fmt = 'png')
pca2d( pca, group= g1 ,show.group.labels = FALSE,axes = "black",show.centroids=TRUE,radius=1.5,show.plane=FALSE,show.ellipses=TRUE,ellipse.ci=0.75,palette=rep(c("steelblue4","aquamarine4","olivedrab3", "orange", "red" )))
rgl.snapshot('2dplot_centroidy.png', fmt = 'png')
#
#  Bez elips ale z centroidami
#
pca3d( pca, group= g1 ,show.group.labels = FALSE,axes = "black",show.centroids=TRUE,radius=1.5,show.plane=FALSE,palette=rep(c("steelblue4","aquamarine4","olivedrab3", "orange", "red" )))
pca2d( pca, group= g1 ,show.group.labels = FALSE,axes = "black",show.centroids=TRUE,radius=1.5,show.plane=FALSE,palette=rep(c("steelblue4","aquamarine4","olivedrab3", "orange", "red" )))
rgl.snapshot('3dplot_centroidy_tylko_bez_elips.png', fmt = 'png')