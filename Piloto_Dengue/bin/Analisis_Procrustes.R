#########################################################################################
###                          Procrustes Analysis                                   ######
#####     Genetic data (E gene and the genome of the Virus Dengue)                 ######
#####     Geographic data (geographic coordinates of countries and localities)     ######
#####                           02- 01 - 18                                        ######
#########################################################################################

# Required Libraries

library(vegan)
library(raster)

# E GENE -- LOCALITY: Data of E gene with reported location

# Genetic distance
Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_loan.csv", header = T, row.names = 1)
Dist_genetica_loan <- as.matrix(Dist_genetica_loan)

# Geographic distance
CoordenadasGE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Coordenadas_loan.csv")
lon_lat <- CoordenadasGE_loan[,c(4,3)]
colnames(lon_lat) <- c("long", "lat")

dist.GE.loan <- pointDistance(lon_lat, lon_lat, lonlat = T, allpairs = T)# distancia geodesica
colnames(dist.GE.loan) <- CoordenadasGE_loan$Location
rownames(dist.GE.loan) <- CoordenadasGE_loan$Location

# Data frame: Accession Number with its respective location
seq_localidad <- data.frame(rownames(Dist_genetica_loan), rownames(dist.GE.loan))

# NMDS for genetics distances E Gene

Mgene <- metaMDS(vegdist(Dist_genetica_loan), autotransform = F, trace = F, trymax = 1000, k=2)
plot(Mgene)

Mgene.stress <- stressplot(Mgene, vegdist(Dist_genetica_loan), main="NMDS Gen E")
ordiplot(Mgene, type = "t")

# NMDS for geographical distances locality

Mgeo <- metaMDS(vegdist(dist.GE.loan), k=2, autotransform = F, trace = F, trymax = 1000, previous.best = dist.GE.loan)
plot(Mgeo)
Mgeo.stress <- stressplot(Mgeo, vegdist(dist.GE.loan), main="NMDS Geografia_Localidades")
ordiplot(Mgeo, type = "t")


#### Procrustes Analysis of NMDS's results
pro.pro <- procrustes(Mgeo,Mgene,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.pro)

plot(pro.pro,kind=2)

####Permutation test of Procrustes analysis

pro.slv4 <- protest(Mgene,Mgeo,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.slv4)

plot(pro.slv4,kind=2,main=NULL)

#---------------------------------------------

# E GEN -- COUNTRY: Data of E gene with reported country

# Genetic distance
Dist_genetica_paan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_paan.csv", header = T, row.names = 1)
Dist_genetica_paan <- as.matrix(Dist_genetica_paan)

# Geographic distance
CoordenadasGE_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Coordenadas_paan.csv")
lon_lat <- CoordenadasGE_paan[,c(3,4)]

dist.GE.paan <- pointDistance(lon_lat, lon_lat, lonlat = T, allpairs = T)
colnames(dist.GE.paan) <- CoordenadasGE_paan$Country
rownames(dist.GE.paan) <- CoordenadasGE_paan$Country

# Data frame: Accession Number with its respective location
seq_pais_E <- data.frame(rownames(Dist_genetica_paan), rownames(dist.GE.paan))

# NMDS for genetics distances E Gene

Mgene <- metaMDS(vegdist(Dist_genetica_paan), autotransform = F, trace = F, trymax = 1000, k=2)
plot(Mgene)

Mgene.stress <- stressplot(Mgene, vegdist(Dist_genetica_paan), main="NMDS Gen E")
ordiplot(Mgene, type = "t")

# NMDS for geographical distances -- country

Mgeo <- metaMDS(vegdist(dist.GE.paan), k=2, autotransform = F, trace = F, trymax = 1000, previous.best = dist.GE.paan)
plot(Mgeo)
Mgeo.stress <- stressplot(Mgeo, vegdist(dist.GE.paan), main="NMDS Geografia_Localidades")
ordiplot(Mgeo, type = "t")

#### Procrustes Analysis of NMDS's results
pro.pro<-procrustes(Mgeo,Mgene,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.pro)

plot(pro.pro,kind=2)

#### Permutation test of Procrustes analysis

pro.slv4 <- protest(Mgene,Mgeo,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.slv4)

plot(pro.slv4,kind=2,main=NULL)

#-----------------------------------------------------------------------------------

# GENOME -- COUNTRY: Data of genome with reported country

# Genetic distance
Dist_genetica_paan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_paan.csv", header = T, row.names = 1)
Dist_genetica_paan <- as.matrix(Dist_genetica_paan)

# Geographic distance
CoordenadasGN_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Coordenadas_paan.csv")
lon_lat <- CoordenadasGN_paan[,c(4,3)]

dist.GN.paan <- pointDistance(lon_lat, lon_lat, lonlat = T, allpairs = T)
colnames(dist.GN.paan) <- CoordenadasGN_paan$Country
rownames(dist.GN.paan) <- CoordenadasGN_paan$Country

# Data frame: Accession Number with its respective location
seq_pais_E <- data.frame(rownames(Dist_genetica_paan), rownames(dist.GN.paan))

# NMDS for genetics distances Genome

Mgene <- metaMDS(vegdist(Dist_genetica_paan), autotransform = F, trace = F, trymax = 1000, k=2)
plot(Mgene)

Mgene.stress <- stressplot(Mgene, vegdist(Dist_genetica_paan), main="NMDS Gen E")
ordiplot(Mgene, type = "t")

# NMDS for geographical distances

Mgeo <- metaMDS(vegdist(dist.GN.paan), k=2, autotransform = F, trace = F, trymax = 1000, previous.best = dist.GN.paan)
plot(Mgeo)
Mgeo.stress <- stressplot(Mgeo, vegdist(dist.GN.paan), main="NMDS Geografia_Localidades")
ordiplot(Mgeo, type = "t")

#### Procrustes Analysis of NMDS's results
pro.pro<-procrustes(Mgeo,Mgene,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.pro)

plot(pro.pro,kind=2)

#### Permutation test of Procrustes analysis

pro.slv4 <- protest(Mgene,Mgeo,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.slv4)

plot(pro.slv4,kind=2,main=NULL)

#-----------------------------------------------------------------------


# GENOME -- LOCALITY: Data of genome with reported locality

# Genetic distance
Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_loan.csv", header = T, row.names = 1)
Dist_genetica_loan <- as.matrix(Dist_genetica_loan)

# Geographic distance
CoordenadasGN_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Coordenadas_loan.csv")
lon_lat <- CoordenadasGN_loan[,c(4,3)]

dist.GN.loan <- pointDistance(lon_lat, lon_lat, lonlat = T, allpairs = T)
colnames(dist.GN.loan) <- CoordenadasGN_loan$Location
rownames(dist.GN.loan) <- CoordenadasGN_loan$Location

# Data frame: Accession Number with its respective location
seq_localidad_E <- data.frame(rownames(Dist_genetica_loan), rownames(dist.GN.loan))

# NMDS for genetics distances Genome

Mgene <- metaMDS(vegdist(Dist_genetica_loan), autotransform = F, trace = F, trymax = 1000, k=2)
plot(Mgene)

Mgene.stress <- stressplot(Mgene, vegdist(Dist_genetica_loan), main="NMDS Gen E")
ordiplot(Mgene, type = "t")

# NMDS for geographical distances

Mgeo <- metaMDS(vegdist(dist.GN.loan), k=2, autotransform = F, trace = F, trymax = 1000, previous.best = dist.GN.loan)
plot(Mgeo)
Mgeo.stress <- stressplot(Mgeo, vegdist(dist.GN.loan), main="NMDS Geografia_Localidades")
ordiplot(Mgeo, type = "t")

#### Procrustes Analysis of NMDS's results
pro.pro<-procrustes(Mgeo,Mgene,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.pro)

plot(pro.pro,kind=2)

#### Permutation test of Procrustes analysis

pro.slv4 <- protest(Mgene,Mgeo,scores="sites",scale=T,symmetric=T,permutation=10000)

plot(pro.slv4)

plot(pro.slv4,kind=2,main=NULL)
