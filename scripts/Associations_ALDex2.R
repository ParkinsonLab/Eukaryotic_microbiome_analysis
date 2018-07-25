####################################################################
#####  ALDex2 - associations between Bacterial OTU and presence
#####  of Eukaryotic microbes or other clinical features



##### ---------- LOAD LIBRARIES, SET DIRECTORIES, IMPORT DATA

library(ALDEx2)

WD <- getwd()
OUTDIR <- paste0(WD, "ALDEx2_analysis")
dir.create(OUTDIR)

DATAFOLDER <- "C:/Users/ProjectData/"


### Import OTU table, taxonomy file and sample metadata

# OTU tables: columns are samples, rows are OTUs
# TAX table: columns are taxonomic ranks, rows are OTUs
# SAMPLEDATA: columns are features of interest, rows are samples/patients
OTU.euk <- read.delim(paste0(DATAFOLDER, "euk_OTU_genus.csv"), sep = "\t", stringsAsFactors = FALSE, row.names = 1)
OTU <- read.delim(paste0(DATAFOLDER, "bact_OTU_genus.txt"), sep = "\t", stringsAsFactors = FALSE, row.names = 1)

TAX <- read.delim(paste0(DATAFOLDER, "OTU_taxonomies.txt"), header=TRUE, row.names = 1, stringsAsFactors = FALSE)

SAMPLEDATA <- read.csv(paste0(DATAFOLDER, "sample_metadata.csv"), row.names = 1, stringsAsFactors = FALSE)
SAMPLEDATA[,"Blastocystis"] <- t(OTU.euk["Blastocystis", rownames(SAMPLEDATA)])
SAMPLEDATA[,"Enterocytozoon"] <- t(OTU.euk["Enterocytozoon", rownames(SAMPLEDATA)])
SAMPLEDATA[,"Cryptosporidium"] <- t(OTU.euk["Cryptosporidium", rownames(SAMPLEDATA)])



##### ---------- CHOOSE/PREPARE VARIABLES AND THRESHOLDS

### Choose variables for testing

var.EUKREADCOUNTS <- c("Blastocystis", "Enterocytozoon", "Cryptosporidium")
var.BINARY <- c("HIV", "Giardia")

SAMPLEDATA.sub <- SAMPLEDATA[,c(var.EUKREADCOUNTS, var.BINARY)]


### Define minimum read cutoffs

# minimum Eukaryote read counts to be considered "positive"
READTHRESHOLD <- c(2,5,50,100)

# minimum reads per bacterial OTU (2 or 5)
READTHRESHOLDBACTERIA <- 2



##### ---------- NORMALIZE BACTERIAL OTU COUNTS (clr) AND ORDINATE

### Choose variables
params <- var.EUKREADCOUNTS
#params <- var.BINARY


### Loop over all thresholds, normalizing data and generating ALDEx2 outputs

# Loop over minimum read thresholds
for(j in 1:length(READTHRESHOLD)) {

  # Loop over features
  for(i in 1:length(params)) {

    # Set conditions and filter OTU table
    PARAMETER <- params[i]
    samples <- rownames(SAMPLEDATA.sub[!is.na(SAMPLEDATA.sub[,PARAMETER]),])
    conds <- SAMPLEDATA.sub[!is.na(SAMPLEDATA.sub[,PARAMETER]),PARAMETER]

    #MINSAMPLETHRESHOLD <- ceiling(0.1*length(samples)) # minimum 10% samples
    MINSAMPLETHRESHOLD <-2 # minimum 2 samples
  
    # Convert feature to binary, if not already
    if(length(params)==length(var.EUKREADCOUNTS)) {
      conds[conds<READTHRESHOLD[j]] <- 0
      conds[conds>=READTHRESHOLD[j]] <- 1
      npos <- length(conds[conds==1])
    }
    
    # Filter out sparse OTUs
    OTU.sub <- OTU[(rowSums(OTU>=READTHRESHOLDBACTERIA))>=MINSAMPLETHRESHOLD, samples]
    
    
    ### ALDEx2 transformation, generates centred log-ratio transformed values

    aldexOUT <- aldex(reads = OTU.sub, conditions = as.character(conds), mc.samples = 128, test = "t", effect = TRUE,
                    include.sample.summary = TRUE, denom = "all", verbose = TRUE)
  
    # Plot
    # we.ep - Expected P-value of Welch's t-test
    png(file = paste0(OUTDIR, "ALDEx2_Bacteria_feature-", PARAMETER, "_subjects", length(samples), "_npos", npos,  "_minread", READTHRESHOLD[j], "_welch.png"), width = 3000, height = 1600, res = 300)
    par(mfrow=c(1,2))
    # MA plot: M (log ratio) vs A (mean average)
    aldex.plot(aldexOUT, type="MA", test="welch") # MA plot
    aldex.plot(aldexOUT, type="MW", test="welch") # MW effect plot
    dev.off()

    # wi.ep - Expected P-value of Wilcoxon rank test
    png(file = paste0(OUTDIR, "ALDEx2_Bacteria_feature-", PARAMETER, "_subjects", length(samples), "_npos", npos,  "_minread", READTHRESHOLD[j], "_wilcox.png"), width = 3000, height = 1600, res = 300)
    par(mfrow=c(1,2))
    # MA plot: M (log ratio) vs A (mean average)
    aldex.plot(aldexOUT, type="MA", test="wilcox", cutoff = 0.05, all.cex = 1, rare.cex = 1, called.cex = 1, called.col = "red") # MA plot
    aldex.plot(aldexOUT, type="MW", test="wilcox", cutoff = 0.05, all.cex = 1, rare.cex = 1, called.cex = 1, called.col = "red") # MW effect plot
    dev.off()
  
    aldexOUT$phylum <- TAX[rownames(aldexOUT),c("phylum")]
    aldexOUT$genus <- TAX[rownames(aldexOUT),c("genus")]
    write.csv(x = aldexOUT, file = paste0(OUTDIR, "ALDEx2_Bacteria_feature-", PARAMETER, "_subjects", length(samples), "_minread", READTHRESHOLD[j], ".csv"), row.names = TRUE)
  }
}



##### ---------- SAVE DATA
save.image(paste0(OUTDIR,"Associations_ALDEx2.rdata"))


# A. Popovic, Parkinson Lab. July 2018