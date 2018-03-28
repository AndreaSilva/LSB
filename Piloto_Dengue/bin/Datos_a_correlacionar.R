###############################################################################################
#####   Selecionar y organizar los datos del Gen E y Genoma que se van a correlacionar   ######
###############################################################################################

##############
##  GEN E  ###
##############

# Datos del Gen E

dat_genE <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/Disimilares97.csv", header = T)

# Secuencias del Gen E

seq_genE <- read.table(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/seq_Disimilaridad97.txt", header= T, sep=" ")

# Datos del gen E que tengan localidad y ano

lo_an <- which(is.na(dat_genE$Location)==F & is.na(dat_genE$Year)==F)

dat_loan <- dat_genE[lo_an,] 

write.csv(dat_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Location_year.csv")

# Secuencias del Gen E que tienen localidad y ano

position_loan <- vector()
for(i in 1:length(dat_loan$N_Accesion)){
  position_loan[i] <- which(seq_genE$id==dat_loan$N_Accesion[i])
}

seq_genE_loan <- seq_genE[position_loan,]

write.table(seq_genE_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Seq_loan.txt")

#-----------------------------------------
##  VECTORES DE COORDENADAS GEOGRAFICAS ##

#GEN E
dge_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Location_year.csv", header = T)
Coordenadas_loan <- dge_loan[,c(12,19,20)]

write.csv(Coordenadas_loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Coordenadas_loan.csv")

#---------------------------------------

## VECTORES DE AÑOS PARA DIST TEMPORAL ##

year.E.loan <- dat_loan$Year
write.csv(year.E.loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_loan.csv")

#----------------------------------

## Gen E con pais y año

# Datos del Gen E

dat_genE <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/Disimilares97.csv", header = T)

# Secuencias del Gen E

seq_genE <- read.table(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/seq_Disimilaridad97.txt", header= T, sep=" ")

# Datos del gen E que tengan pais y ano

lo_an <- which(is.na(dat_genE$Country)==F & is.na(dat_genE$Year)==F)

dat_loan <- dat_genE[lo_an,] 

write.csv(dat_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Country_year.csv")

# Secuencias del Gen E que tienen pais y ano

position_loan <- vector()
for(i in 1:length(dat_loan$N_Accesion)){
  position_loan[i] <- which(seq_genE$id==dat_loan$N_Accesion[i])
}

seq_genE_loan <- seq_genE[position_loan,]

write.table(seq_genE_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Seq_paan.txt")

#-----------------------------------------
##  VECTORES DE COORDENADAS GEOGRAFICAS ##

#GEN E
dge_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Country_year.csv", header = T)
Coordenadas_loan <- dge_loan[,c(11,16,15)]

write.csv(Coordenadas_loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Coordenadas_paan.csv")

#---------------------------------------

## VECTORES DE AÑOS PARA DIST TEMPORAL ##

year.E.loan <- dat_loan$Year
write.csv(year.E.loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_paan.csv")

###################
###  GENOMA    ####
###################

## 1. Genoma del virus del Dengue con Localidad y Año

# Datos de Genoma
dat_genoma <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/Disimilares97.csv", header = T)

# Secuencias del Genoma
seq_genoma <- read.table(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/Seq_disimilaridad97.txt", header= T, sep=" ")

# Datos de Genoma que tiene localidad y ano
Genomalo_an <- which(is.na(dat_genoma$Location)==F & is.na(dat_genoma$Year)==F & is.na(dat_genoma$lat_location)==F)
dat_Genomalo_an <- dat_genoma[Genomalo_an,]
dat_Genomalo_an <- dat_Genomalo_an[-25,]

write.csv(dat_Genomalo_an, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Location_year.csv")

# Secuencias de Genoma que tiene localidad y ano
posinoma_loan <- vector()
for(i in 1:length(dat_Genomalo_an$N_Accesion)){
  posinoma_loan[i] <- which(seq_genoma$id==dat_Genomalo_an$N_Accesion[i])
}

seq_genoma_loan <- seq_genoma[posinoma_loan,]

write.table(seq_genoma_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Seq_loan.txt")

#-----------------------------

## VECTORESDE COORDENAS GEOGRAFICAS

#GENOMA
# LOCALIDAD Y ANO
CoordenadasG_loan <- dat_Genomalo_an[,c(11,18,19)]
write.csv(CoordenadasG_loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Coordenadas_loan.csv")

#------------------------------

## VECTORES DE AÑO 

year_loan <- dat_Genomalo_an$Year
write.csv(year_loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_loan.csv")


## 2. Genoma del virus del Dengue con Pais y Año

# Datos del genoma que tienen pais y ano
Genomapa_an <- which(is.na(dat_genoma$Country)==F & is.na(dat_genoma$Year)==F)
dat_Genoma_paan <- dat_genoma[Genomapa_an,]
write.csv(dat_Genoma_paan, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Country_year.csv")

# Secuencias del genoma con pais y ano
posinoma_paan <- vector()
for(i in 1:length(dat_Genoma_paan$N_Accesion)){
  posinoma_paan[i] <- which(seq_genoma$id==dat_Genoma_paan$N_Accesion[i])
}

seq_genoma_paan <- seq_genoma[posinoma_paan,]

write.table(seq_genoma_paan, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Seq_pa_an.txt")

#--------------------------
## Vectores de coordenadas geograficas ##

#PAIS Y ANO 
CoordenadasG_paan <- dat_Genoma_paan[,c(10,14,15)]
write.csv(CoordenadasG_paan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Coordenadas_paan.csv")

#----------------------------------

## VECTORES DE AÑO ##

year_paan <- dat_Genoma_paan$Year
write.csv(year_paan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_paan.csv")
