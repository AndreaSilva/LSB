

disttemporal <- function(vector){
  distT <- matrix(0,length(vector),length(vector))
  for (i in 1:length(vector)){
    for (j in 1:length(vector)){
      distT[i,j] <- abs(vector[i]-vector[j])
    }
  }
  return(distT)
}

# GEN E #
# Gen E localidad y ano

yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_loan.csv")
yearE_loan <- yearE_loan$x

distE_loan <- disttemporal(yearE_loan)
colnames(distE_loan) <- yearE_loan
rownames(distE_loan) <- yearE_loan
write.csv(distE_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_temporalea/yearE_loan.csv")

# Gen E dengue 1 
yearD1_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/yearD1_loan.csv")
yearD1_loan <- yearD1_loan$x

distD1_loan <- disttemporal(yearD1_loan)
colnames(distD1_loan) <- yearD1_loan
rownames(distD1_loan) <- yearD1_loan
write.csv(distD1_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_temporalea/yearD1_loan.csv")

# Gen E dengue 2 
yearD2_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/yearD2_loan.csv")
yearD2_loan <- yearD2_loan$x

distD2_loan <- disttemporal(yearD2_loan)
colnames(distD2_loan) <- yearD2_loan
rownames(distD2_loan) <- yearD2_loan
write.csv(distD2_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_temporalea/yearD2_loan.csv")

# Gen E dengue 3 
yearD3_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/yearD3_loan.csv")
yearD3_loan <- yearD3_loan$x

distD3_loan <- disttemporal(yearD3_loan)
colnames(distD3_loan) <- yearD3_loan
rownames(distD3_loan) <- yearD3_loan
write.csv(distD3_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_temporalea/yearD3_loan.csv")

# Gen E dengue 4 
yearD4_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/yearD4_loan.csv")
yearD4_loan <- yearD4_loan$x

distD4_loan <- disttemporal(yearD4_loan)
colnames(distD4_loan) <- yearD4_loan
rownames(distD4_loan) <- yearD4_loan
write.csv(distD4_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_temporalea/yearD4_loan.csv")
#-----------------------------

# Gen E pais y ano

yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_paan.csv")
yearE_loan <- yearE_loan$x

distE_loan <- disttemporal(yearE_loan)
colnames(distE_loan) <- yearE_loan
rownames(distE_loan) <- yearE_loan
write.csv(distE_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_temporalea/yearE_paan.csv")

# Gen E dengue 1 
yearD1_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/yearD1_paan.csv")
yearD1_loan <- yearD1_loan$x

distD1_loan <- disttemporal(yearD1_loan)
colnames(distD1_loan) <- yearD1_loan
rownames(distD1_loan) <- yearD1_loan
write.csv(distD1_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_temporalea/yearD1_paan.csv")

# Gen E dengue 2 
yearD2_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/yearD2_paan.csv")
yearD2_loan <- yearD2_loan$x

distD2_loan <- disttemporal(yearD2_loan)
colnames(distD2_loan) <- yearD2_loan
rownames(distD2_loan) <- yearD2_loan
write.csv(distD2_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_temporalea/yearD2_paan.csv")

# Gen E dengue 3 
yearD3_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/yearD3_paan.csv")
yearD3_loan <- yearD3_loan$x

distD3_loan <- disttemporal(yearD3_loan)
colnames(distD3_loan) <- yearD3_loan
rownames(distD3_loan) <- yearD3_loan
write.csv(distD3_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_temporalea/yearD3_paan.csv")

# Gen E dengue 4 
yearD4_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/yearD4_paan.csv")
yearD4_loan <- yearD4_loan$x

distD4_loan <- disttemporal(yearD4_loan)
colnames(distD4_loan) <- yearD4_loan
rownames(distD4_loan) <- yearD4_loan
write.csv(distD4_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_temporalea/yearD4_paan.csv")

#----------------------------------------------------------------------------------------------------
# GENOMA #
# Genoma localidad y ano

yearG_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_loan.csv")
yearG_loan <- yearG_loan$x

distG_loan <- disttemporal(yearG_loan)
colnames(distG_loan) <- yearG_loan
rownames(distG_loan) <- yearG_loan
write.csv(distG_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_temporales/distG_loan.csv")

# Genoma Dengue 1 localidad y ano

yearD1_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/yearD1_loan.csv")
yearD1_loan <- yearD1_loan$x

distD1_loan <- disttemporal(yearD1_loan)
colnames(distD1_loan) <- yearD1_loan
rownames(distD1_loan) <- yearD1_loan
write.csv(distD1_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_temporales/distD1_loan.csv")

# Genoma Dengue 2 localidad y ano

yearD2_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/yearD2_loan.csv")
yearD2_loan <- yearD2_loan$x

distD2_loan <- disttemporal(yearD2_loan)
colnames(distD2_loan) <- yearD2_loan
rownames(distD2_loan) <- yearD2_loan
write.csv(distD2_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_temporales/distD2_loan.csv")

# Genoma Dengue 3 localidad y ano

yearD3_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/yearD3_loan.csv")
yearD3_loan <- yearD3_loan$x

distD3_loan <- disttemporal(yearD3_loan)
colnames(distD3_loan) <- yearD3_loan
rownames(distD3_loan) <- yearD3_loan
write.csv(distD3_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_temporales/distD3_loan.csv")

# Genoma Dengue 4 localidad y ano

yearD4_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/yearD4_loan.csv")
yearD4_loan <- yearD4_loan$x

distD4_loan <- disttemporal(yearD4_loan)
colnames(distD4_loan) <- yearD4_loan
rownames(distD4_loan) <- yearD4_loan
write.csv(distD4_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_temporales/distD4_loan.csv")
#-------------------

# Genoma pais y ano

yearG_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_paan.csv")
yearG_paan <- yearG_paan$x

distG_paan <- disttemporal(yearG_paan)
colnames(distG_paan) <- yearG_paan
rownames(distG_paan) <- yearG_paan
write.csv(distG_paan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_temporales/distG_paan.csv")

# Genoma Dengue 1 pais y ano

yearD1_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/yearD1_paan.csv")
yearD1_paan <- yearD1_paan$x

distD1_paan <- disttemporal(yearD1_paan)
colnames(distD1_paan) <- yearD1_paan
rownames(distD1_paan) <- yearD1_paan
write.csv(distD1_paan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_temporales/distD1_paan.csv")

# Genoma Dengue 2 pais y ano

yearD2_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/yearD2_paan.csv")
yearD2_paan <- yearD2_paan$x

distD2_paan <- disttemporal(yearD2_paan)
colnames(distD2_paan) <- yearD2_paan
rownames(distD2_paan) <- yearD2_paan
write.csv(distD2_paan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_temporales/distD2_paan.csv")

# Genoma Dengue 3 pais y ano

yearD3_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/yearD3_paan.csv")
yearD3_paan <- yearD3_paan$x

distD3_paan <- disttemporal(yearD3_paan)
colnames(distD3_paan) <- yearD3_paan
rownames(distD3_paan) <- yearD3_paan
write.csv(distD3_paan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_temporales/distD3_paan.csv")

# Genoma Dengue 4 pais y ano

yearD4_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/yearD4_paan.csv")
yearD4_paan <- yearD4_paan$x

distD4_paan <- disttemporal(yearD4_paan)
colnames(distD4_paan) <- yearD4_paan
rownames(distD4_paan) <- yearD4_paan
write.csv(distD4_paan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_temporales/distD4_paan.csv")

