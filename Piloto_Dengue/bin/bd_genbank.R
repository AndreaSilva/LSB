#######################################################################
####           Obtain Dengue database from GenBank                #####
####         For the Dengue virus, using the rentrez package      #####
####                                 2016                         #####
#######################################################################


#libraries required
library(rentrez)

dengue <- "Dengue virus[organism]"

# search the "dengue" records in the NCBI nucleotide database
dengue_search <- entrez_search(db="nuccore", term=dengue, retmax=18389 ,use_history = TRUE)

#general information of the search
dengue_search

#Total number of records found
dengue_search$count

#ids of the searched sequences
dengue_search$ids
dengue_search$ids[1:30]

#see what each search contains
dengue_summs_1 <- entrez_summary(db="nuccore", id=dengue_search$ids)

# vector with additional information that I want to obtain from each sequence
ater <- c("title", "organism", "extra", "caption", "slen", "subname")

#caption: Numeros de accesion
#organism: organismo tipo de dengue
#title: titlo del gb
#extra: da el numero de gi y gb
#taxid: xref= taxon
#crearedate: fecha de creacion del genbank
#slen: longitud de la secuencia
#subtype: "genotipo, pais, fuente de aislado, fecha de colecta
#subname: describe a subtype

#caption: Accession numbers
#organism: type organism of dengue
#title: titlo del gb
#extra: gives the number of gi and gb
#taxid: xref = taxon
#crearedate: date of creation of the genbank
#slen: length of the sequence
#subtype: "genotype, country, source of isolated, date of collection
#subname: describes subtype

# obtener informacion adicional que describí en el vector anterior
ejm_1 <- extract_from_esummary(dengue_summs_1,ater, simplify = T)
unname(ejm)

# Loop to put all the information in a dataframe
datos <- data.frame()

for(i in dengue_search$ids[16760:18389]){
  dengue_summs <- entrez_summary(db="nuccore", id=i)
  ejm <- extract_from_esummary(dengue_summs, ater, simplify = T)
  datos[i,1] <- unname(ejm)[[1]]
  datos[i,2] <- unname(ejm)[[2]]
  datos[i,3] <- unname(ejm)[[3]]
  datos[i,4] <- unname(ejm)[[4]]
  datos[i,5] <- unname(ejm)[[5]]
  datos[i,6] <- unname(ejm)[[6]]
  #Sys.sleep()
}

# Save downloaded data
write.csv(datos, file="/home/andrea/PROYECTO_LSB/Piloto_Dengue/data/compendio_datos.csv")
write.table(datos, file="/home/andrea/PROYECTO_LSB/Piloto_Dengue/data/total_datos.csv", sep = ";")

###########################################################################################
# NO LO USE
#Paquete seqinr
#library(seqinr)

#Elegir la base de datos en la que quiero buscar
#choosebank("genbank")

#buscar por nombre de la especie"SP"
#seqdengue <- query("seqdengue","SP=Dengue virus")

#atributos de la busqueda
#attributes(seqdengue)

#numero de datos encontrados
#seqdengue$nelem

#informacion de los datos encontrados: numero de accsesion
#datos <- seqdengue$req

#cerrar conexion
#closebank()

###Guardar datos (N° accesion y longitud) de las secuencias

#crear dataframe vacio
#datos2 <- data.frame()

#bucle para llenar el dataframe con datos de las secuencias
#for(i in 1:18424){
#  datos2[i,"accesion"] <- datos[[i]][1]
#  datos2[i, "length"] <- attr(datos[[i]],"length")
#}

#guardar datos en un .csv
#write.table(datos2, file="seqdengue.csv", sep=",")

#annots <- getAnnot(datos[[1]])

#annots[1:20]=="ORGANISM "
