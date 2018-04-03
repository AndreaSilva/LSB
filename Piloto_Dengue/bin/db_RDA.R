#####################################################################################
####           Distance-based redundancy analysis (db-RDA)                       ####
###          For genetic data (gene E and genome), geography                     ####
###  (lat-long of countries and localities) and temporary (years) of Dengue virus####
######################################################################################

# library required

library(vegan)

##############
## E GENE ####
##############

# 1. E GenE - locality - year

# genetic matrix
Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_loan.csv", header = T, row.names = 1)
Dist_genetica_loan <- as.matrix(Dist_genetica_loan)
Dist_genetica_loan <- Dist_genetica_loan+0.15 # eliminar valores imaginarios
for(i in 1:nrow(Dist_genetica_loan)){
  Dist_genetica_loan[i,i] <- 0
}

# explanatory variables:

#years
yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_loan.csv")
yearE_loan <- yearE_loan[2:2]
colnames(yearE_loan) <-"Years"

# geography
CoordenadasGE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Coordenadas_loan.csv")
lon_lat <- CoordenadasGE_loan[,c(4,3)]

geotem <- cbind(lon_lat,yearE_loan)

# db-RDA
dbrda_temporal <- capscale(formula = Dist_genetica_loan ~. , data = geotem)

plot(dbrda_temporal, main="")
ordiplot(dbrda_temporal,type="t",scaling=2)

dbrda_temporal$CCA$tot.chi/dbrda_temporal$tot.chi

# permutation test

anova.cca(dbrda_temporal)
anova(dbrda_temporal)
anova(dbrda_temporal, by="terms", perm.max=999)
anova(dbrda_temporal, by="axis", perm=999)
anova(dbrda_temporal, by="margin", perm=999)

#--------------------------------------

#2. Gen E - paises - years

# genetic matrix
Dist_genetica_paan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_paan.csv", header = T, row.names = 1)
Dist_genetica_paan <- as.matrix(Dist_genetica_paan)
Dist_genetica_paan <- Dist_genetica_paan+0.65
for(i in 1:nrow(Dist_genetica_paan)){
  Dist_genetica_paan[i,i] <- 0
}

# explanatory variables:

#years

yearE_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_paan.csv")
yearE_paan <- yearE_paan[2:2]
colnames(yearE_paan) <-"Years"

#geography
CoordenadasGE_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Coordenadas_paan.csv")
lon_lat <- CoordenadasGE_paan[,c(3,4)]

geotem <- cbind(lon_lat, yearE_paan)

# db-RDA
dbrda_temporal <- capscale(formula = Dist_genetica_paan ~. , data = geotem)

plot(dbrda_temporal, type="t")
ordiplot(dbrda_temporal,type="t",scaling=2)

dbrda_temporal$CCA$tot.chi/dbrda_temporal$tot.chi

# permutation test

anova.cca(dbrda_temporal)
anova(dbrda_temporal)
anova(dbrda_temporal, by="terms", perm.max=999)
anova(dbrda_temporal, by="axis", perm=999)

###-----------------------------------------------------

##############
### GENOME ###
##############

# 3. Genome - locality - years

# genetic matrix
Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_loan.csv", header = T, row.names = 1)
Dist_genetica_loan <- as.matrix(Dist_genetica_loan)
Dist_genetica_loan <- Dist_genetica_loan+0.01
for(i in 1:nrow(Dist_genetica_loan)){
  Dist_genetica_loan[i,i] <- 0
}

# explanatory variables:
#years
yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_loan.csv")
yearE_loan <- yearE_loan[2:2]
colnames(yearE_loan) <-"Years"

#geography
CoordenadasGN_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Coordenadas_loan.csv")
lon_lat <- CoordenadasGN_loan[,c(4,3)]

geotem <- cbind(lon_lat, yearE_loan)

# db-RDA
dbrda_temporal <- capscale(formula = Dist_genetica_loan ~. , data = geotem)

plot(dbrda_temporal, main="")
ordiplot(dbrda_temporal,type="t",scaling=2)

dbrda_temporal$CCA$tot.chi/dbrda_temporal$tot.chi

# permutation test

anova.cca(dbrda_temporal)
anova(dbrda_temporal)
anova(dbrda_temporal, by="terms", perm.max=999)
anova(dbrda_temporal, by="axis", perm=999)
#-----------------------------------------

#4. Genome - country - years

# genetic matrix
Dist_genetica_paan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_paan.csv", header = T, row.names = 1)
Dist_genetica_paan <- as.matrix(Dist_genetica_paan)
Dist_genetica_paan <- Dist_genetica_paan+0.1
for(i in 1:nrow(Dist_genetica_paan)){
  Dist_genetica_paan[i,i] <- 0
}

# explanatory variables:

#years
yearE_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_paan.csv")
yearE_paan <- yearE_paan[2:2]
colnames(yearE_paan) <-"Years"

#geography
CoordenadasGN_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Coordenadas_paan.csv")
lon_lat <- CoordenadasGN_paan[,c(4,3)]

geotem <- cbind(lon_lat,yearE_paan)

# db-RDA
dbrda_temporal <- capscale(formula = Dist_genetica_paan ~. , data = geotem)

plot(dbrda_temporal, type="t")
ordiplot(dbrda_temporal,type="t",scaling=2)

dbrda_temporal$CCA$tot.chi/dbrda_temporal$tot.chi

# permutation test

anova.cca(dbrda_temporal)
anova(dbrda_temporal)
anova(dbrda_temporal, by="terms", perm.max=999)
anova(dbrda_temporal, by="axis", perm=999)

#---------------------------------------------------------------------

#library(vegan)

# Gen E

# genetic matrix
#Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_loan.csv", header = T, row.names = 1)
#Dist_genetica_loan <- as.matrix(Dist_genetica_loan)
#Dist_genetica_loan <- Dist_genetica_loan+0.2
#for(i in 1:nrow(Dist_genetica_loan)){
#  Dist_genetica_loan[i,i] <- 0
#}

# explanatory variables: years
#yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_loan.csv")
#yearE_loan <- yearE_loan[2:2]
#colnames(yearE_loan) <-"Years"

# db-RDA
#dbrda_temporal <- capscale(formula = Dist_genetica_loan ~. , data = yearE_loan)

#plot(dbrda_temporal, main="")

#--------------------------------------

# Gen E

# genetic matrix
#Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_paan.csv", header = T, row.names = 1)
#Dist_genetica_loan <- as.matrix(Dist_genetica_loan)
#Dist_genetica_loan <- Dist_genetica_loan+0.65
#for(i in 1:nrow(Dist_genetica_loan)){
#  Dist_genetica_loan[i,i] <- 0
#}

# explanatory variables: years
#yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_paan.csv")
#yearE_loan <- yearE_loan[2:2]
#colnames(yearE_loan) <-"Years"

# db-RDA
#dbrda_temporal <- capscale(formula = Dist_genetica_loan ~. , data = yearE_loan)

#plot(dbrda_temporal, type="t")

###-----------------------------------------------------

# Genoma

# genetic matrix
#Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_loan.csv", header = T, row.names = 1)
#Dist_genetica_loan <- as.matrix(Dist_genetica_loan)
#Dist_genetica_loan <- Dist_genetica_loan+0.01
#for(i in 1:nrow(Dist_genetica_loan)){
#  Dist_genetica_loan[i,i] <- 0
#}

# explanatory variables: years
#yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_loan.csv")
#yearE_loan <- yearE_loan[2:2]
#colnames(yearE_loan) <-"Years"

# db-RDA
#dbrda_temporal <- capscale(formula = Dist_genetica_loan ~. , data = yearE_loan)

#plot(dbrda_temporal, main="")


# Genoma pais

# genetic matrix
#Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_paan.csv", header = T, row.names = 1)
#Dist_genetica_loan <- as.matrix(Dist_genetica_loan)
#Dist_genetica_loan <- Dist_genetica_loan+0.1
#for(i in 1:nrow(Dist_genetica_loan)){
#  Dist_genetica_loan[i,i] <- 0
#}

# explanatory variables: years
#yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_paan.csv")
#yearE_loan <- yearE_loan[2:2]
#colnames(yearE_loan) <-"Years"

# db-RDA
#dbrda_temporal <- capscale(formula = Dist_genetica_loan ~. , data = yearE_loan)

#plot(dbrda_temporal, type="t")
