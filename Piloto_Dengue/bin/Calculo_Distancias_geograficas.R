#############################################################################
####            Euclidean geographical distance calculation           #######
#####                  Raster -- point distances                      #######
#############################################################################

# Required library
library(raster)

## E GEN ##

# Gen E -- locality: Data of E gene with reported location
CoordenadasGE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Coordenadas_loan.csv")
lon_lat <- CoordenadasGE_loan[,c(4,3)]

dist.GE.loan <- pointDistance(lon_lat, lon_lat, lonlat = F, allpairs = T)
colnames(dist.GE.loan) <- CoordenadasGE_loan$Location
rownames(dist.GE.loan) <- CoordenadasGE_loan$Location
write.csv(dist.GE.loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_geograficas/dist.GE.loan.csv")

#--------------------------------

# E Gen -- Country: Data of E gene with reported country
CoordenadasGE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Coordenadas_paan.csv")
lon_lat <- CoordenadasGE_loan[,c(3,4)]

dist.GE.loan <- pointDistance(lon_lat, lon_lat, lonlat = F, allpairs = T)
colnames(dist.GE.loan) <- CoordenadasGE_loan$Country
rownames(dist.GE.loan) <- CoordenadasGE_loan$Country
write.csv(dist.GE.loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_geograficas/dist.GE.paan.csv")

## GENOME ##

# Genome -- locality: Data of genome with reported location
CoordenadasGN_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Coordenadas_loan.csv")
lon_lat <- CoordenadasGN_loan[,c(4,3)]

dist.GN.loan <- pointDistance(lon_lat, lon_lat, lonlat = F, allpairs = T)
colnames(dist.GN.loan) <- CoordenadasGN_loan$Location
rownames(dist.GN.loan) <- CoordenadasGN_loan$Location
write.csv(dist.GN.loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_geograficas/dist.GN.loan.csv")

#-----------------------

# Genome -- country: Data of genome with reported country

CoordenadasGN_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Coordenadas_paan.csv")
lon_lat <- CoordenadasGN_paan[,c(4,3)]

dist.GN.paan <- pointDistance(lon_lat, lon_lat, lonlat = F, allpairs = T)
colnames(dist.GN.paan) <- CoordenadasGN_paan$Country
rownames(dist.GN.paan) <- CoordenadasGN_paan$Country
write.csv(dist.GN.paan,file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_geograficas/dist.GN.paan.csv")
