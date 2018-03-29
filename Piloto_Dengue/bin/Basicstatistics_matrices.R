#######################################################################################
###                Basic statistics of genetic distance matrices              #########
###                            (maximum, minimum and mean)                    #########
#######################################################################################

# Matrix E gene with locality
mgenEl <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_loan.csv", header = T, row.names = 1)
mgenEl <- as.matrix(mgenEl)
ugenEl <- mgenEl[-which(upper.tri(mgenEl)==F)]
max(ugenEl)
min(ugenEl)
mean(ugenEl)

#Matrix E gene with country
mgenEp <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genE/Fraccional_paan.csv", header = T, row.names = 1)
mgenEp <- as.matrix(mgenEp)
ugenEp <- mgenEp[-which(upper.tri(mgenEp)==F)]
max(ugenEp)
min(ugenEp)
mean(ugenEp)

#Matrix genome with locality
mgenOl <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_loan.csv", header = T, row.names = 1)
mgenOl <- as.matrix(mgenOl)
ugenOl <- mgenOl[-which(upper.tri(mgenOl)==F)]
max(ugenOl)
min(ugenOl)
mean(ugenOl)

#Matrix genome with country
mgenOp <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Matrices_distancias/matrices_genoma/Fraccional_paan.csv", header = T, row.names = 1)
mgenOp <- as.matrix(mgenOp)
ugenOp <- mgenOp[-which(upper.tri(mgenOp)==F)]
max(ugenOp)
min(ugenOp)
mean(ugenOp)

