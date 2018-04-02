###############################################################################################
#####   Selecionar y organizar los datos del Gen E y Genoma que se van a correlacionar   ######
###############################################################################################

##############
##  E GENE  ###
##############

# E Gene Data

dat_genE <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/Disimilares97.csv", header = T)

# E Gene sequences

seq_genE <- read.table(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/seq_Disimilaridad97.txt", header= T, sep=" ")

# Data of E gene with location and year

lo_an <- which(is.na(dat_genE$Location)==F & is.na(dat_genE$Year)==F)

dat_loan <- dat_genE[lo_an,]

write.csv(dat_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Location_year.csv")

# E gene sequences with location and year

position_loan <- vector()
for(i in 1:length(dat_loan$N_Accesion)){
  position_loan[i] <- which(seq_genE$id==dat_loan$N_Accesion[i])
}

seq_genE_loan <- seq_genE[position_loan,]

write.table(seq_genE_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Seq_loan.txt")

#-----------------------------------------
##  GEOGRAPHICAL COORDINATE VECTORS ##

# E GENE
dge_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Location_year.csv", header = T)
Coordenadas_loan <- dge_loan[,c(12,19,20)]

write.csv(Coordenadas_loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Coordenadas_loan.csv")

#---------------------------------------

## YEARS VECTORS - TEMPORARY DISTANCE ##

year.E.loan <- dat_loan$Year
write.csv(year.E.loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_loan.csv")

#----------------------------------

## E Gene with country and year

# E gene Data

dat_genE <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/Disimilares97.csv", header = T)

# E gene sequences

seq_genE <- read.table(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/seq_Disimilaridad97.txt", header= T, sep=" ")

# Datos del gen E que tengan pais y ano

lo_an <- which(is.na(dat_genE$Country)==F & is.na(dat_genE$Year)==F)

dat_loan <- dat_genE[lo_an,]

write.csv(dat_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Country_year.csv")

# E gene sequences with country and year

position_loan <- vector()
for(i in 1:length(dat_loan$N_Accesion)){
  position_loan[i] <- which(seq_genE$id==dat_loan$N_Accesion[i])
}

seq_genE_loan <- seq_genE[position_loan,]

write.table(seq_genE_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Seq_paan.txt")

#-----------------------------------------
##  GEOGRAPHICAL COORDINATE VECTORS ##

#GEN E
dge_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Country_year.csv", header = T)
Coordenadas_loan <- dge_loan[,c(11,16,15)]

write.csv(Coordenadas_loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Coordenadas_paan.csv")

#---------------------------------------

## YEARS VECTORS - TEMPORARY DISTANCE ##

year.E.loan <- dat_loan$Year
write.csv(year.E.loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_paan.csv")

###################
###  GENOME    ####
###################

## Genome Data with locality and year

# Genoma Data
dat_genoma <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/Disimilares97.csv", header = T)

# Genome sequences
seq_genoma <- read.table(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/Seq_disimilaridad97.txt", header= T, sep=" ")

# Genoma Data with locality and year
Genomalo_an <- which(is.na(dat_genoma$Location)==F & is.na(dat_genoma$Year)==F & is.na(dat_genoma$lat_location)==F)
dat_Genomalo_an <- dat_genoma[Genomalo_an,]
dat_Genomalo_an <- dat_Genomalo_an[-25,]

write.csv(dat_Genomalo_an, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Location_year.csv")

# Genome sequences with locality and year
posinoma_loan <- vector()
for(i in 1:length(dat_Genomalo_an$N_Accesion)){
  posinoma_loan[i] <- which(seq_genoma$id==dat_Genomalo_an$N_Accesion[i])
}

seq_genoma_loan <- seq_genoma[posinoma_loan,]

write.table(seq_genoma_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Seq_loan.txt")

#-----------------------------

## GEOGRAPHICAL COORDINATE VECTORS

#GENOME
# LOCALITY AND YEAR
CoordenadasG_loan <- dat_Genomalo_an[,c(11,18,19)]
write.csv(CoordenadasG_loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Coordenadas_loan.csv")

#------------------------------

## YEARS VECTORS - TEMPORARY DISTANCE

year_loan <- dat_Genomalo_an$Year
write.csv(year_loan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_loan.csv")

# Genome Data with country and year
Genomapa_an <- which(is.na(dat_genoma$Country)==F & is.na(dat_genoma$Year)==F)
dat_Genoma_paan <- dat_genoma[Genomapa_an,]
write.csv(dat_Genoma_paan, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Country_year.csv")

# Genome sequences with country and year
posinoma_paan <- vector()
for(i in 1:length(dat_Genoma_paan$N_Accesion)){
  posinoma_paan[i] <- which(seq_genoma$id==dat_Genoma_paan$N_Accesion[i])
}

seq_genoma_paan <- seq_genoma[posinoma_paan,]

write.table(seq_genoma_paan, file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Seq_pa_an.txt")

#--------------------------
## GEOGRAPHICAL COORDINATE VECTORS ##

#country and year 
CoordenadasG_paan <- dat_Genoma_paan[,c(10,14,15)]
write.csv(CoordenadasG_paan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Coordenadas_paan.csv")

#----------------------------------

## YEARS VECTORS - TEMPORARY DISTANCE

year_paan <- dat_Genoma_paan$Year
write.csv(year_paan,file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_paan.csv")
