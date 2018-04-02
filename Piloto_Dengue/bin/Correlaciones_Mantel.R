########################################################################################
#####                                 Mantel test                                   ####
####  genetic data (E gene and genome) and geographic data (countries and localities) ##
#####                    and temporary (years) of the Dengue virus                  ####
########################################################################################

## E GENE Vs GEOGRAPHIC Vs TIME (LOCALITY AND YEARS)

# kmer distances
Fraccional_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_loan.csv", header = T, row.names = 1)

# Geographic distance
distGE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_geograficas/dist_GE_loan.csv")
distGE_loan$X <- colnames(distGE_loan[2:30])
rownames(distGE_loan) <- distGE_loan$X
distGE_loan <- distGE_loan[2:30]

#Tempoaral distance
yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_temporalea/yearE_loan.csv")
yearE_loan$X <- colnames(yearE_loan[2:30])
rownames(yearE_loan) <- yearE_loan$X
yearE_loan <- yearE_loan[2:30]

####################
## Mantel Test###
####################

# Required library
library(vegan)

mantel(xdis= Fraccional_loan, ydis= distGE_loan)
mantel(xdis= Fraccional_loan, ydis= yearE_loan)
mantel(xdis= distGE_loan, ydis= yearE_loan)

#-----------------------------------------------------------------------------------------------------

## E GENE Vs GEOGRAPHIC Vs TIME (COUNTRY AND YEARS)

# Kmer distance
Fraccional_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_paan.csv", header = T, row.names = 1)

# Geographic distance
distGE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_geograficas/dist.GE.paan.csv")
distGE_loan$X <- colnames(distGE_loan[2:131])
rownames(distGE_loan) <- distGE_loan$X
distGE_loan <- distGE_loan[2:131]

#Tempoaral distance
yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/matrices_temporalea/yearE_paan.csv")
yearE_loan$X <- colnames(yearE_loan[2:131])
rownames(yearE_loan) <- yearE_loan$X
yearE_loan <- yearE_loan[2:131]

####################
## Mantel Test###
####################

# Required library
library(vegan)

mantel(xdis= Fraccional_loan, ydis= distGE_loan)
mantel(xdis= Fraccional_loan, ydis= yearE_loan)
mantel(xdis= distGE_loan, ydis= yearE_loan)

#------------------------------------------------------------------------------------------------------
## GENOME Vs GEOGRAPHIC Vs TIME (LOCALITY AND YEARS)

# Kmer distance
Fraccional_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_loan.csv", header = T, row.names = 1)

# Geographic distance
distG_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_geograficas/dist.GN.loan.csv")
distG_loan$X <- colnames(distG_loan[2:42])
rownames(distG_loan) <- distG_loan$X
distG_loan <- distG_loan[2:42]

#Tempoaral distance
yearG_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_temporales/distG_loan.csv")
yearG_loan$X <- colnames(yearG_loan[2:42])
rownames(yearG_loan) <- yearG_loan$X
yearG_loan <- yearG_loan[2:42]

####################
## Mantel Test###
####################

library(vegan)

mantel(xdis= Fraccional_loan, ydis= distG_loan)
mantel(xdis= Fraccional_loan, ydis= yearG_loan)
mantel(xdis= distG_loan, ydis= yearG_loan)

#-----------------------------

## GENOME Vs GEOGRAPHIC Vs TIME (COUNTRY AND YEARS)

# Kmer distance
Fraccional_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_paan.csv", header = T, row.names = 1)

# Geographic distance
distG_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_geograficas/dist.GN.paan.csv")
distG_paan$X <- colnames(distG_paan[2:84])
rownames(distG_paan) <- distG_paan$X
distG_paan <- distG_paan[2:84]

#Tempoaral distance
yearG_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/matrices_temporales/distG_paan.csv")
yearG_paan$X <- colnames(yearG_paan[2:84])
rownames(yearG_paan) <- yearG_paan$X
yearG_paan <- yearG_paan[2:84]

##################
## Mantel Test ###
##################

library(vegan)

mantel(xdis= Fraccional_paan, ydis= distG_paan)
mantel(xdis= Fraccional_paan, ydis= yearG_paan)
mantel(xdis= distG_paan, ydis= yearG_paan)
