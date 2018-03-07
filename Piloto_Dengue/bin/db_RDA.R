##########################################################################################
####         Analisis de redundancia canonica basado en distancias (db-RDA)           #### 
###              Para datos genéticos (gen E y genoma), geograficos                   ####
###     (lat-long de paises y localidades) y temporales (años) del virus del Dengue   ####
##########################################################################################

# libreria requerida

library(vegan)

##############
## GEN E  ####
##############

# 1. Gen E - localidades - años

# matriz genetica
Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_loan.csv", header = T, row.names = 1)
Dist_genetica_loan <- as.matrix(Dist_genetica_loan)
Dist_genetica_loan <- Dist_genetica_loan+0.15 # eliminar valores imaginarios
for(i in 1:nrow(Dist_genetica_loan)){
  Dist_genetica_loan[i,i] <- 0
}

# variables explicatorias: 

#años
yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_loan.csv")
yearE_loan <- yearE_loan[2:2]
colnames(yearE_loan) <-"Years"

# geografia
CoordenadasGE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Coordenadas_loan.csv")
lon_lat <- CoordenadasGE_loan[,c(4,3)]

geotem <- cbind(lon_lat,yearE_loan)

# db-RDA
dbrda_temporal <- capscale(formula = Dist_genetica_loan ~. , data = geotem)

plot(dbrda_temporal, main="")
ordiplot(dbrda_temporal,type="t",scaling=2)

dbrda_temporal$CCA$tot.chi/dbrda_temporal$tot.chi

# test de permutacion

anova.cca(dbrda_temporal)
anova(dbrda_temporal)
anova(dbrda_temporal, by="terms", perm.max=999)
anova(dbrda_temporal, by="axis", perm=999)
anova(dbrda_temporal, by="margin", perm=999)

#--------------------------------------

#2. Gen E - paises - años

# matriz genetica
Dist_genetica_paan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_paan.csv", header = T, row.names = 1)
Dist_genetica_paan <- as.matrix(Dist_genetica_paan)
Dist_genetica_paan <- Dist_genetica_paan+0.65
for(i in 1:nrow(Dist_genetica_paan)){
  Dist_genetica_paan[i,i] <- 0
}

# variables explicatorias: 

#años

yearE_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_paan.csv")
yearE_paan <- yearE_paan[2:2]
colnames(yearE_paan) <-"Years"

#geografia
CoordenadasGE_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/Coordenadas_paan.csv")
lon_lat <- CoordenadasGE_paan[,c(3,4)]

geotem <- cbind(lon_lat, yearE_paan)

# db-RDA
dbrda_temporal <- capscale(formula = Dist_genetica_paan ~. , data = geotem)

plot(dbrda_temporal, type="t")
ordiplot(dbrda_temporal,type="t",scaling=2)

dbrda_temporal$CCA$tot.chi/dbrda_temporal$tot.chi

# test de permutacion

anova.cca(dbrda_temporal)
anova(dbrda_temporal)
anova(dbrda_temporal, by="terms", perm.max=999)
anova(dbrda_temporal, by="axis", perm=999)

###-----------------------------------------------------

##############
### GENOMA ###
##############

# 3. Genoma - localidad - años

# matriz genetica
Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_loan.csv", header = T, row.names = 1)
Dist_genetica_loan <- as.matrix(Dist_genetica_loan)
Dist_genetica_loan <- Dist_genetica_loan+0.01
for(i in 1:nrow(Dist_genetica_loan)){
  Dist_genetica_loan[i,i] <- 0
}

# variables explicatorias: 
#años
yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_loan.csv")
yearE_loan <- yearE_loan[2:2]
colnames(yearE_loan) <-"Years"

#geografia
CoordenadasGN_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Coordenadas_loan.csv")
lon_lat <- CoordenadasGN_loan[,c(4,3)]

geotem <- cbind(lon_lat, yearE_loan)

# db-RDA
dbrda_temporal <- capscale(formula = Dist_genetica_loan ~. , data = geotem)

plot(dbrda_temporal, main="")
ordiplot(dbrda_temporal,type="t",scaling=2)

dbrda_temporal$CCA$tot.chi/dbrda_temporal$tot.chi

# test de permutacion

anova.cca(dbrda_temporal)
anova(dbrda_temporal)
anova(dbrda_temporal, by="terms", perm.max=999)
anova(dbrda_temporal, by="axis", perm=999)
#-----------------------------------------

#4. Genoma - pais - años

# matriz genetica
Dist_genetica_paan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_paan.csv", header = T, row.names = 1)
Dist_genetica_paan <- as.matrix(Dist_genetica_paan)
Dist_genetica_paan <- Dist_genetica_paan+0.1
for(i in 1:nrow(Dist_genetica_paan)){
  Dist_genetica_paan[i,i] <- 0
}

# variables explicatorias: 

#años
yearE_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_paan.csv")
yearE_paan <- yearE_paan[2:2]
colnames(yearE_paan) <-"Years"

#geografia
CoordenadasGN_paan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/Coordenadas_paan.csv")
lon_lat <- CoordenadasGN_paan[,c(4,3)]

geotem <- cbind(lon_lat,yearE_paan)

# db-RDA
dbrda_temporal <- capscale(formula = Dist_genetica_paan ~. , data = geotem)

plot(dbrda_temporal, type="t")
ordiplot(dbrda_temporal,type="t",scaling=2)

dbrda_temporal$CCA$tot.chi/dbrda_temporal$tot.chi

# test de permutacion

anova.cca(dbrda_temporal)
anova(dbrda_temporal)
anova(dbrda_temporal, by="terms", perm.max=999)
anova(dbrda_temporal, by="axis", perm=999)

#---------------------------------------------------------------------

#library(vegan)

# Gen E

# matriz genetica
#Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_loan.csv", header = T, row.names = 1)
#Dist_genetica_loan <- as.matrix(Dist_genetica_loan)
#Dist_genetica_loan <- Dist_genetica_loan+0.2
#for(i in 1:nrow(Dist_genetica_loan)){
#  Dist_genetica_loan[i,i] <- 0
#}

# variables explicatorias: años
#yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_loan.csv")
#yearE_loan <- yearE_loan[2:2]
#colnames(yearE_loan) <-"Years"

# db-RDA
#dbrda_temporal <- capscale(formula = Dist_genetica_loan ~. , data = yearE_loan)

#plot(dbrda_temporal, main="")

#--------------------------------------

# Gen E

# matriz genetica
#Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_paan.csv", header = T, row.names = 1)
#Dist_genetica_loan <- as.matrix(Dist_genetica_loan)
#Dist_genetica_loan <- Dist_genetica_loan+0.65
#for(i in 1:nrow(Dist_genetica_loan)){
#  Dist_genetica_loan[i,i] <- 0
#}

# variables explicatorias: años
#yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/sub_database/year_paan.csv")
#yearE_loan <- yearE_loan[2:2]
#colnames(yearE_loan) <-"Years"

# db-RDA
#dbrda_temporal <- capscale(formula = Dist_genetica_loan ~. , data = yearE_loan)

#plot(dbrda_temporal, type="t")

###-----------------------------------------------------

# Genoma 

# matriz genetica
#Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_loan.csv", header = T, row.names = 1)
#Dist_genetica_loan <- as.matrix(Dist_genetica_loan)
#Dist_genetica_loan <- Dist_genetica_loan+0.01
#for(i in 1:nrow(Dist_genetica_loan)){
#  Dist_genetica_loan[i,i] <- 0
#}

# variables explicatorias: años
#yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_loan.csv")
#yearE_loan <- yearE_loan[2:2]
#colnames(yearE_loan) <-"Years"

# db-RDA
#dbrda_temporal <- capscale(formula = Dist_genetica_loan ~. , data = yearE_loan)

#plot(dbrda_temporal, main="")


# Genoma pais

# matriz genetica
#Dist_genetica_loan <- read.csv("/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_paan.csv", header = T, row.names = 1)
#Dist_genetica_loan <- as.matrix(Dist_genetica_loan)
#Dist_genetica_loan <- Dist_genetica_loan+0.1
#for(i in 1:nrow(Dist_genetica_loan)){
#  Dist_genetica_loan[i,i] <- 0
#}

# variables explicatorias: años
#yearE_loan <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genoma_completo/sub_database/year_paan.csv")
#yearE_loan <- yearE_loan[2:2]
#colnames(yearE_loan) <-"Years"

# db-RDA
#dbrda_temporal <- capscale(formula = Dist_genetica_loan ~. , data = yearE_loan)

#plot(dbrda_temporal, type="t")