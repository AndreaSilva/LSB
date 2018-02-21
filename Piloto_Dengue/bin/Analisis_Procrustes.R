#########################################################################################
###            Análisis de procrustes para datos genéticos (gen E y genoma)        ######
#####     y geográficos (coordenadas de paises y localidades)del virus del Dengue   #####
###                                 02- 01 - 18                                     #####
#########################################################################################

# librerias requeridas

library(vegan)
library(raster)

# GEN E LOCALIDAD

# 1. gen E localidad

# distancia genetica
Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_loan.csv", header = T, row.names = 1)
Dist_genetica_loan <- as.matrix(Dist_genetica_loan)

# distancia geografica
CoordenadasGE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Coordenadas_loan.csv")
lon_lat <- CoordenadasGE_loan[,c(4,3)]
colnames(lon_lat) <- c("long", "lat")

dist.GE.loan <- pointDistance(lon_lat, lon_lat, lonlat = T, allpairs = T)
colnames(dist.GE.loan) <- CoordenadasGE_loan$Location
rownames(dist.GE.loan) <- CoordenadasGE_loan$Location

# Data frame de numero de acceso con su respectivo ubicacion
seq_localidad <- data.frame(rownames(Dist_genetica_loan), rownames(dist.GE.loan))

# NMDS de distancias geneticas Gen E

Mgene <- metaMDS(vegdist(Dist_genetica_loan), autotransform = F, trace = F, trymax = 1000, k=2)
plot(Mgene)

Mgene.stress <- stressplot(Mgene, vegdist(Dist_genetica_loan), main="NMDS Gen E")
ordiplot(Mgene, type = "t")

# NMDS de distancias geograficas localidades gen e

Mgeo <- metaMDS(vegdist(dist.GE.loan), k=2, autotransform = F, trace = F, trymax = 1000, previous.best = dist.GE.loan)
plot(Mgeo)
Mgeo.stress <- stressplot(Mgeo, vegdist(dist.GE.loan), main="NMDS Geografia_Localidades")
ordiplot(Mgeo, type = "t")


#### Procrustes
pro.pro <- procrustes(Mgeo,Mgene,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.pro)

plot(pro.pro,kind=2)

####Permutation test of Procrustes analysis

pro.slv4 <- protest(Mgene,Mgeo,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.slv4)

plot(pro.slv4,kind=2,main=NULL)

#---------------------------------------------

# GEN E PAIS

# 1. gen E pais

# distancia genetica
Dist_genetica_paan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_paan.csv", header = T, row.names = 1)
Dist_genetica_paan <- as.matrix(Dist_genetica_paan)

# distancia geografica
CoordenadasGE_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Coordenadas_paan.csv")
lon_lat <- CoordenadasGE_paan[,c(3,4)]

dist.GE.paan <- pointDistance(lon_lat, lon_lat, lonlat = T, allpairs = T)
colnames(dist.GE.paan) <- CoordenadasGE_paan$Country
rownames(dist.GE.paan) <- CoordenadasGE_paan$Country

# Data frame de numero de acceso con su respectivo ubicacion
seq_pais_E <- data.frame(rownames(Dist_genetica_paan), rownames(dist.GE.paan))

# NMDS de distancias geneticas Gen E

Mgene <- metaMDS(vegdist(Dist_genetica_paan), autotransform = F, trace = F, trymax = 1000, k=2)
plot(Mgene)

Mgene.stress <- stressplot(Mgene, vegdist(Dist_genetica_paan), main="NMDS Gen E")
ordiplot(Mgene, type = "t")

# NMDS de distancias geograficas paises gen e

Mgeo <- metaMDS(vegdist(dist.GE.paan), k=2, autotransform = F, trace = F, trymax = 1000, previous.best = dist.GE.paan)
plot(Mgeo)
Mgeo.stress <- stressplot(Mgeo, vegdist(dist.GE.paan), main="NMDS Geografia_Localidades")
ordiplot(Mgeo, type = "t")

#### Procrustes
pro.pro<-procrustes(Mgeo,Mgene,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.pro)

plot(pro.pro,kind=2)

####Permutation test of Procrustes analysis

pro.slv4 <- protest(Mgene,Mgeo,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.slv4)

plot(pro.slv4,kind=2,main=NULL)

#-----------------------------------------------------------------------------------

# GENOMA PAIS

# 1. genoma pais

# distancia genetica
Dist_genetica_paan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_paan.csv", header = T, row.names = 1)
Dist_genetica_paan <- as.matrix(Dist_genetica_paan)

# distancia geografica
CoordenadasGN_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Coordenadas_paan.csv")
lon_lat <- CoordenadasGN_paan[,c(4,3)]

dist.GN.paan <- pointDistance(lon_lat, lon_lat, lonlat = T, allpairs = T)
colnames(dist.GN.paan) <- CoordenadasGN_paan$Country
rownames(dist.GN.paan) <- CoordenadasGN_paan$Country

# Data frame de numero de acceso con su respectivo ubicacion
seq_pais_E <- data.frame(rownames(Dist_genetica_paan), rownames(dist.GN.paan))

# NMDS de distancias geneticas Gen E

Mgene <- metaMDS(vegdist(Dist_genetica_paan), autotransform = F, trace = F, trymax = 1000, k=2)
plot(Mgene)

Mgene.stress <- stressplot(Mgene, vegdist(Dist_genetica_paan), main="NMDS Gen E")
ordiplot(Mgene, type = "t")

# NMDS de distancias geograficas paises gen e

Mgeo <- metaMDS(vegdist(dist.GN.paan), k=2, autotransform = F, trace = F, trymax = 1000, previous.best = dist.GN.paan)
plot(Mgeo)
Mgeo.stress <- stressplot(Mgeo, vegdist(dist.GN.paan), main="NMDS Geografia_Localidades")
ordiplot(Mgeo, type = "t")

#### Procrustes
pro.pro<-procrustes(Mgeo,Mgene,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.pro)

plot(pro.pro,kind=2)

####Permutation test of Procrustes analysis

pro.slv4 <- protest(Mgene,Mgeo,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.slv4)

plot(pro.slv4,kind=2,main=NULL)

#-----------------------------------------------------------------------


# GENOMA LOCALIDAD

# 1. genoma localidad

# distancia genetica
Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_loan.csv", header = T, row.names = 1)
Dist_genetica_loan <- as.matrix(Dist_genetica_loan)

# distancia geografica
CoordenadasGN_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Coordenadas_loan.csv")
lon_lat <- CoordenadasGN_loan[,c(4,3)]

dist.GN.loan <- pointDistance(lon_lat, lon_lat, lonlat = T, allpairs = T)
colnames(dist.GN.loan) <- CoordenadasGN_loan$Location
rownames(dist.GN.loan) <- CoordenadasGN_loan$Location

# Data frame de numero de acceso con su respectivo ubicacion
seq_localidad_E <- data.frame(rownames(Dist_genetica_loan), rownames(dist.GN.loan))

# NMDS de distancias geneticas Gen E

Mgene <- metaMDS(vegdist(Dist_genetica_loan), autotransform = F, trace = F, trymax = 1000, k=2)
plot(Mgene)

Mgene.stress <- stressplot(Mgene, vegdist(Dist_genetica_loan), main="NMDS Gen E")
ordiplot(Mgene, type = "t")

# NMDS de distancias geograficas paises gen e

Mgeo <- metaMDS(vegdist(dist.GN.loan), k=2, autotransform = F, trace = F, trymax = 1000, previous.best = dist.GN.loan)
plot(Mgeo)
Mgeo.stress <- stressplot(Mgeo, vegdist(dist.GN.loan), main="NMDS Geografia_Localidades")
ordiplot(Mgeo, type = "t")

#### Procrustes
pro.pro<-procrustes(Mgeo,Mgene,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.pro)

plot(pro.pro,kind=2)

####Permutation test of Procrustes analysis

pro.slv4 <- protest(Mgene,Mgeo,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.slv4)

plot(pro.slv4,kind=2,main=NULL)
