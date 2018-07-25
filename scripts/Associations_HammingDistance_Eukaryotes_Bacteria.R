############################################################################################
#####     Hamming distance - Eukaryotes (18S V4+V5 or transITS amplicons) vs. Bacteria (16S)



##### ---------- LOAD LIBRARIES, SET DIRECTORIES, IMPORT DATA


### Load libraries
library("e1071")
library("gplots")
library("reshape2")


### Set directories
WD <- getwd()
OUTDIR <- paste0(WD, "Hamming_BetweenBlastocystisGroups_18S/")
dir.create(OUTDIR)

FILEDIRECTORY <- "C:/Users/ProjectData/"


### Import OTU and taxonomy tables (at different taxonomic levels)

TAXLEVELS <- c("genus", "class", "phylum")
BACTAXLEVEL <- TAXLEVELS[1]
EUKTAXLEVEL <- TAXLEVELS[1]

Bac_OTU <- read.delim(paste0(FILEDIRECTORY, "16S/bact_OTU_", BACTAXLEVEL, ".txt"), sep = "\t", stringsAsFactors = FALSE, row.names = 1)
Bac_TAX <- read.delim(paste0(FILEDIRECTORY, "16S/OTU_taxonomies.txt"), header=TRUE, stringsAsFactors = FALSE, row.names = 1)
Bac_OTU$genus <- TAX[rownames(Bac_OTU),"genus"]

v.Euk_OTU <- read.csv(file = paste0(FILEDIRECTORY, "18S/Euk_OTU_", EUKTAXLEVEL, ".csv"),  header = TRUE, stringsAsFactors = FALSE, row.names = 1)
v.Euk_TAX <- read.delim(paste0(FILEDIRECTORY, "18S/OTU_taxonomies.txt"), header=TRUE, stringsAsFactors = FALSE, row.names = 1)
v.Euk_OTU$genus <- v.Euk_TAX[rownames(v.Euk_OTU),"genus"]

t.Euk_OTU <- read.delim(file = paste0(FILEDIRECTORY, "transITS/Euk_OTU_", EUKTAX, ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
t.Euk_TAX <- read.delim(paste0(FILEDIRECTORY, "transITS/OTU_taxonomies.txt"), header=TRUE, stringsAsFactors = FALSE, row.names = 1)
t.Euk_OTU$genus <- t.Euk_TAX[rownames(t.Euk_OTU),"genus"]



##### ---------- SET PARAMETERS, PREPARE/FILTER OTU TABLES

# Choose OTU tables
BACOTU <- Bac_OTU
EUKOTU <- v.Euk_OTU

# Replace OTU number with taxonomic name
rownames(BACOTU) <- Bac_TAX[rownames(BACOTU), BACTAXLEVEL]
rownames(EUKOTU) <- Euk_TAX[rownames(EUKOTU), EUKTAXLEVEL]

# Remove samples not found in Eukaryotic data, ensure sample orders
BACOTU <- BACOTU[,colnames(BACOTU) %in% colnames(EUKOTU)]
BACOTU <- BACOTU[,colnames(EUKOTU)]

# Subset to OTUs present in 20%-80%, 15%-85% or 10%-90% of samples, with X minimum reads
MIN <- 0.2
MAX <- 0.8
MINREAD <- c(2,5,20,50,100)



##### ---------- CALCULATE HAMMING DISTANCES


### Loop over minimum read read thresholds
for(i in 1:length(MINREAD)) {

  
  # Prefilter OTU tables
  B <- BACOTU[(rowMeans(BACOTU >= MINREAD[i])>=MIN)&(rowMeans(BACOTU >= MINREAD[i])<=MAX),]
  E <- EUKOTU[(rowMeans(EUKOTU >= MINREAD[i])>=MIN)&(rowMeans(EUKOTU >= MINREAD[i])<=MAX),]


  # Combine matrices and convert to binary
  allotu <- rbind(E, B)
  allotu[allotu < MINREAD[i]] <- 0
  allotu[allotu >= MINREAD[i]] <- 1

  
  ### Calculate Hamming distances
  otu_hamdist <- hamming.distance(as.matrix(allotu))
  otu_hamdist <- as.data.frame(otu_hamdist)

  
  ### Plot heatmap of Hamming distances
  my_palette <- colorRampPalette(c("red", "white", "blue"))
  heatmap.2(as.matrix(otu_hamdist), margins = c(13, 13), cexRow=1.1, cexCol = 1.1, tracecol=NA, col = my_palette, dendrogram='none', key=FALSE, lwid=c(0.1,4), lhei=c(0.1,4))

  # Subset heatmap to only show eukaryotes in rows, bacteria in columns
  sub <- otu_hamdist[rownames(E),rownames(B)]
  heatmap.2(as.matrix(sub), margins = c(11, 13), cexCol = 1.1, cexRow=1.1, tracecol=NA, col = my_palette, dendrogram='none', key=FALSE, lwid=c(0.1,4), lhei=c(0.1,4))
  
  
  ### Write Hamming distance output
  write.csv(sub, file = paste0(OUTDIR, "HammingDistanceTable_", Sys.Date(), "_Euk", EUKTAXLEVEL, "_Bac", BACTAXLEVEL, "_", MIN, "-", MAX, ".csv"), row.names = TRUE)

}



##### ---------- SAVE DATA

save.image(paste0(OUTDIR,"Hamming_Euk_Bacteria_associations.rdata"))


# A. Popovic, Parkinson Lab. July 2018