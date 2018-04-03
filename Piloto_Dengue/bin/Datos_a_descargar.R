################################################################################
####                      Download complete Genome                       #######
###              and of the E gene sequences of the Dengue virus         #######
###                        from the Accession Number                     #######
################################################################################

## Download complete Genome Sequences ###

# Read file with the database
bd_dengue <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/data/Base_Datos/Base_Datos_Dengue/bd_Dengue.csv", stringsAsFactors = F)

# Data positions with the following characteristics (to eliminate them):

# 1. Data with NA in the "country" column and the "year" column (In both at the same time)
country_year_na <- which(is.na(bd_dengue$Country) & is.na(bd_dengue$Year))

# 2. Data with unknown serotype
dsserotype <- which(bd_dengue$Serotype=="DENV")

# 3. Data that are clones
dsclone <- which(bd_dengue$Serotype=="Clone")

# 4. Unverified data
dsunverified <- which(bd_dengue$Serotype=="UNVERIFIED")

# 5. Data that are chimeras
dschimeric <- which(bd_dengue$Serotype=="Chimeric")

# Join data positions to eliminate
bd_eliminar <- c(country_year_na,dsserotype,dschimeric, dsclone,dsunverified)

bd_eliminar <- unique(bd_eliminar)

# Delete data that is not downloaded
datos <- bd_dengue[-bd_eliminar,]

# New data base with E gene and genome data
write.csv(datos,file = "/home/andrea/LSB/Piloto_Dengue/data/Base_Datos/Base_Datos_Dengue/bd_Dengue_v02.csv")

# From the new database: Extract data from the complete genome
Complete_genome <- datos[which(datos$Gene=="Complete_Genome" & datos$Size_sequence>10100),]

# List of Accession Aumbers
ID <- Complete_genome$N_Accesion

which(ID=="AJ487271")
No.Accesion <-ID[3893:3916]
ID[3893]
length(ID)
----------------------------------------------------------------------------------------------------
# These sequences were downloaded using the function written by Viviana Romero-Alarcon
# in the code downloadCDSgb.R
# executed with the code CDS_dowland_genbank.R
# https://github.com/alarconvv/download.CDS.GenBank


#Some sequences were not downloaded for the following reasons:

# The CDS is not reported in GenBank
# An access number is not reported, only the gb number is reported and this number is not recognized by the code
# It has a length less than 1100 bp
# They are labeled by GenBank as a clone, chimera, unverified
# Does not report the serotype

# only the gb number is reported:
#  "NC_002640" "NC_001474" "NC_001477" "NC_001475"

# The CDS is not reported in GenBank
#sincds <- c("EU920839", "KC131141", "KC131140", "HM631854", "HM631853" ,"HM631852", "HM631851", "HM631855",
#            "FJ913015", "JQ045685", "JQ045683", "JQ045682", "JQ045681", "JQ045680", "JQ045679", "JQ045678",
#            "JQ045677", "JQ045676", "JQ045675", "JQ045674", "JQ045673", "JQ045672", "JQ045671", "JQ045669",
#            "JQ045627", "M87512", "EU179861", "EU179860", "FJ177308")

# Clone
# FJ828987, FJ828986, AY145122, AY145121, AY145123, U88537, U88536, U88535, AB543624
# AF375822, AF326827, AF326826, AF326825, DQ285561,
----------------------------------------------------------------------------------------------------
# Partition to run by parts.

ID <- ID[1:5]
ID <- ID[501:1000]
ID <- ID[1001:1500]
ID <- ID[1501:2000]
ID <- ID[2001:2500]
ID <- ID[2501:3000]
ID <- ID[3001:3500]
ID <- ID[3501: 4000]
ID <- ID[4001:4016]

sequences_dengue <- read.csv(file = "/home/andrea/LSB/Piloto_Dengue/bin/descargar CDS r/Dengue_total")

#-----------------------------------------------------------------------------------------

###################################
# Download E gene sequences  #
####################################

# Required libraries

library(ape)
library(seqinr)

# Extract the E gene data with a sequence length between
# 1300 to 1600 bp, taking into account that the length of the E gene reported by the
# GenBank reference sequences are 1485 bp

Gen_E <- datos[datos$Gene=="E" & datos$Size_sequence>1300 & datos$Size_sequence<1600, ]

# List of Accession Numbers

ids <- Gen_E$N_Accesion
length(ids)
which(ids=="AB111090")

# The search of the sequences from the access number

genE_sequences <- read.GenBank(ids)

genE_sequences_ids <- paste(attr(genE_sequences, "species"), names(genE_sequences), sep = "_")

# Save the sequences

write.dna(genE_sequences,"/home/andrea/LSB/Piloto_Dengue/data/Secuencias_descargadas/Secuencias_genE/Gen_E.fasta", format = "fasta")
