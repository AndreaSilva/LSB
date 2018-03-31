########################################################################################
####          Calculation of the absolute delta between years                     ######
####              data: years of Dengue virus registration                         #####
####                                   15 - 12 - 17                                #####
########################################################################################

# Function for the calculation of absolute delta between years

disttemporal <- function(vector){
  distT <- matrix(0,length(vector),length(vector))
  for (i in 1:length(vector)){
    for (j in 1:length(vector)){
      distT[i,j] <- abs(vector[i]-vector[j])
    }
  }
  return(distT)
}

# Data: years of data registration of E Gene and Genome of Dengue virus
# E GEN #
# E Gen whit locality and year

yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_loan.csv")
yearE_loan <- yearE_loan$x

distE_loan <- disttemporal(yearE_loan)
colnames(distE_loan) <- yearE_loan
rownames(distE_loan) <- yearE_loan
write.csv(distE_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_temporalea/yearE_loan.csv")

#-----------------------------

# E Gen country and year

yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_paan.csv")
yearE_loan <- yearE_loan$x

distE_loan <- disttemporal(yearE_loan)
colnames(distE_loan) <- yearE_loan
rownames(distE_loan) <- yearE_loan
write.csv(distE_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_temporalea/yearE_paan.csv")

#----------------------------------------------------------------------------------------------------
# GENOME #
# Genome whit locality and year

yearG_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_loan.csv")
yearG_loan <- yearG_loan$x

distG_loan <- disttemporal(yearG_loan)
colnames(distG_loan) <- yearG_loan
rownames(distG_loan) <- yearG_loan
write.csv(distG_loan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_temporales/distG_loan.csv")

#-------------------

# Genome whit country and year

yearG_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_paan.csv")
yearG_paan <- yearG_paan$x

distG_paan <- disttemporal(yearG_paan)
colnames(distG_paan) <- yearG_paan
rownames(distG_paan) <- yearG_paan
write.csv(distG_paan, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_temporales/distG_paan.csv")
