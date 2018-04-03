#######################################################################
##     Conteo de k-mer y calculo de las distancias geneticas  #########
####     distancia Euclidiana, de Mahalanobis y Conteo comun  #########
#####                       de k-mer                          #########
#######################################################################

# Data reading: Dengue DNA sequences

#-----------------------------------------------------------------------------------
## E GENE ##

#1. Sequences of E Gene with location and year
dat_seq <- read.table(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Seq_loan.txt", header= T, sep=" ")

seq <- dat_seq$seq

#-------------------

#1. Sequences of E Gene with country and year
dat_seq <- read.table(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Seq_paan.txt", header= T, sep=" ")

seq <- dat_seq$seq

#------------------------------------------------------------------------------------------

## GENOME ##

#1. Sequences of Genome with location and year
dat_seq <- read.table(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Seq_loan.txt", header= T, sep=" ")
seq <- dat_seq$seq

#-------

#2. Sequences of Genome with country and year
dat_seq <- read.table(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Seq_pa_an.txt", header= T, sep=" ")
seq <- dat_seq$seq

# Count of k-mers for triplets k = 3

library(stringi)

kmer3 <-  paste(rep(c("A","C","G","T"),each=16),rep(c("A","C","G","T"),each=4),c("A","C","G","T"),sep = "")

length(kmer3)

result <- t(sapply(seq, stri_count_fixed, pattern=kmer3, overlap=TRUE))
colnames(result) <- kmer3
result <- as.data.frame(result)

######################################################################################################
########################
# Euclidian distance #
# sum(p(s1)-p(s2))^2  ##
########################

Euclidian <- function(table){
  distE <- matrix(0,nrow(table),nrow(table))
  for (i in 1:nrow(table)){
    for (j in 1:nrow(table)){
      distE[i,j] <- sum((table[i,]-table[j,])^2)
    }
  }
  return(distE)
}

Euclidiana <- Euclidian(result)

#-----------------------------------------------------------------------------------------------
## E GENE MATRIX ##

# matrix genE con localidad y ano
colnames(Euclidiana) <- dat_seq$id
rownames(Euclidiana) <- dat_seq$id
Dist_Euclidiana <- write.csv(Euclidiana, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Euclidiana_loan.csv")

#------------

# matrix genE con pais y ano
colnames(Euclidiana) <- dat_seq$id
rownames(Euclidiana) <- dat_seq$id
Dist_Euclidiana <- write.csv(Euclidiana, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Euclidiana_paan.csv")

#---------------------------------------------------------------------------------------------

## GENOME MATRIX

# matrix genoma con localidad y ano
colnames(Euclidiana) <- dat_seq$id
rownames(Euclidiana) <- dat_seq$id
Dist_Euclidiana <- write.csv(Euclidiana, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Euclidiana_loan.csv")

#----------

# matrix genoma con pais y ano
colnames(Euclidiana) <- dat_seq$id
rownames(Euclidiana) <- dat_seq$id
Dist_Euclidiana <- write.csv(Euclidiana, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Euclidiana_paan.csv")

#######################################################################################################
###################
# Dist Mahalanobis #
####################

#Calculo de capacidad de solapamiento para todas las tripletas posibles
#Calculation of overlapping capacity for all possible triples

library(reshape2)
kmerspt <- colsplit(kmer3, "", names=c(1,2,3))

overlap <- data.frame()
for(i in 1:nrow(kmerspt)){
  for(j in 1:ncol(kmerspt)){
    overlap[i,1] <- as.numeric(kmerspt[i,3]==kmerspt[i,1])
    overlap[i,2] <- as.numeric(kmerspt[i,3]==kmerspt[i,2])
  }
}

rownames(overlap) <- kmer3

#Calculo de la longitud de cada una de las secuencias
# N ---> las longitudes de todas las secuencias
#Calculation of the length of each of the sequences

lonN <- data.frame()

for (i in 1:length(seq)){
  lonN[i,1] <- length(strsplit(as.character(seq[i]), NULL)[[1]])
}


# Funcion de la varianza para k = 3
# Var[f(a1..ak)] = E (1-1/4^k) - 2/4^2k (k-1)(N-(3/2K)+1) + 2/4^k sum(N-K+1-t) Jt/4^t
# Function of the variance for k = 3

Varian <- function(table1, table2, k=3){
  vari <- data.frame()
  for(i in 1:nrow(table1)){
    for(j in 1:nrow(table2)){
      vari[i,j] <- (((table1[i,1]-k+1)*(1/4)^k) * (1-((table1[i,1]-k+1)*(1/4)^k)) + ((((1/4)^k)^2)*((table1[i,1]-k+1)-k)*((table1[i,1]-k+1)-k+1)) + (2*(1/4)^k) * sum(((table1[i,1]-k+1)-1)*(sum(table2[j,]))*((1/4)^1)))
    }
  }
  return(vari)
}

varianza <- Varian(lonN,overlap)
colnames(varianza) <- kmer3

# Calculo de la Distancia de mahalanobis
# Mahalanobis distance

nobis <- data.frame()
for (i in 1:nrow(result)){
  for (j in 1:nrow(varianza)){
    nobis[i,j] <- sum(((result[i,]/varianza[i,])-(result[j,]/varianza[j,]))^2)
  }
}

Mahalanobis <- nobis

#--------------------------------------------------------------------------------------------

## MATRICES GEN E

# Matrix de Gen E con localidad y ano
colnames(Mahalanobis) <- dat_seq$id
rownames(Mahalanobis) <- dat_seq$id
Dist_Mahalanobis <- write.csv(Mahalanobis, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Mahalanobis_loan.csv")
#-----------------------------------

# Matrix de Gen E con pais y ano
colnames(Mahalanobis) <- dat_seq$id
rownames(Mahalanobis) <- dat_seq$id
Dist_Mahalanobis <- write.csv(Mahalanobis, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Mahalanobis_paan.csv")

#--------------------------------------------------------------------------------------------------

## MATRICES GENOMA

# matrix genoma con localidad y ano
colnames(Mahalanobis) <- dat_seq$id
rownames(Mahalanobis) <- dat_seq$id
Dist_Mahalanobis <- write.csv(Mahalanobis, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Mahalanobis_loan.csv")

#----------------
# matrix genoma con pais y ano
colnames(Mahalanobis) <- dat_seq$id
rownames(Mahalanobis) <- dat_seq$id
Dist_Mahalanobis <- write.csv(Mahalanobis, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Mahalanobis_paan.csv")

# 3:25.82 min

########################################################################################################
#########################################################
# Fractional Common k-mer Count                      ##
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC373290/  ##
# 1 - (sum (min(p(s1),p(s2))/(min(l1,l2)-k +1)))       ##
# min(p(s1),p(s2) --> valor minimo del conteo de k-mers #
# min(l1,l2) --> Valor minimo de longitud de secuencia ##
# k --> longitud de subsecuencia                       ##
# log --> transformacion logaritmica                   ##
#########################################################

# Calculo de los minimos entre las frecuencias por pares de secuencias
# Calculation of the minimums between the frequencies by pairs of sequences

soul_sal <- data.frame()

for(i in 1:ncol(result)){
  soul <- vector()
  for(f in 1:length(result[,i])){
    for(ty in 1:length(result[,i])){
      soul <- c(soul, min(c(result[f,i], result[ty,i])))
    }
  }
  soul_sal[1:441,i] <- soul
}

sume <- t(apply(soul_sal, 1, sum))
sume <- as.data.frame(t(sume))
##40 seg
# Calculo de lo minimos entre las longitudes por pares de secuencias
# Calculation of the minimum between the lengths by pairs of sequences

minlon <- vector()
for (i in 1:nrow(lonN)){
  algo <- 1:length(lonN[,1])
  for(n in algo){
    minlon <- c(minlon,min(lonN[i,1],lonN[n,1]))
  }
}
min_lon <- as.data.frame(minlon)

# Calculo de la ecuacion completa
# Fraccional k-mer common

Fractional <- function(table, table1, k=3){ # sume, minlon
  commonk <- data.frame()
  for(j in 1:nrow(table)){
    commonk[1,j] <- 1- (table[j,1]/(table1[j,1]-k+1))
  }
  return(commonk)
}

Fractional_common <- Fractional(sume, min_lon)

# 15 seg

#Construcción de la matriz para Fraccional k-mer common
# Construction of the matrix for Fractional k-mer common

mmtx <- matrix(data = Fractional_common, nrow = 21, ncol = 21, byrow = FALSE,
       dimnames = NULL)

#-------------------------------------------------------------------------------------------------

## MATRICES GEN E

# Matrix gen E con localidad y ano
colnames(mmtx) <- dat_seq$id
rownames(mmtx) <- dat_seq$id
Dist_Fractional <- write.csv(mmtx, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_loan.csv")

#-------------------

# Matrix gen E con pais y ano
colnames(mmtx) <- dat_seq$id
rownames(mmtx) <- dat_seq$id
Dist_Fractional <- write.csv(mmtx, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_paan.csv")

#------------------------------------------------------------------------------------------------

## MATRICES GENOMA ##

# Matrix genoma con localidad y ano
colnames(mmtx) <- dat_seq$id
rownames(mmtx) <- dat_seq$id
Dist_Fractional <- write.csv(mmtx, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_loan.csv")
#-----

# Matrix genoma con pais y ano
colnames(mmtx) <- dat_seq$id
rownames(mmtx) <- dat_seq$id
Dist_Fractional <- write.csv(mmtx, file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_paan.csv")
