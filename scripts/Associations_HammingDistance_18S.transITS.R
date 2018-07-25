#######################################################################
#####     Hamming distance - 18S V4+V5 and transITS amplicons     #####



##### ---------- LOAD LIBRARIES, SET DIRECTORIES, IMPORT DATA

library("e1071")
library("gplots")
library("reshape2")
library("ggplot2")


### Set directories

WD <- getwd()
OUTDIR <- paste0(WD, "18S-transITS_HammingDistance/")
dir.create(OUTDIR)

FILEDIRECTORY <- "C:/Users/ProjectDirectory/"


### Import OTU tables

# 18S data
v.phylum <- read.csv(file = paste0(FILEDIRECTORY, "18SV4V5/Euk_OTU_phylum.csv"), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
v.genus <- read.csv(file = paste0(FILEDIRECTORY, "18SV4V5/Euk_OTU_genus.csv"),  header = TRUE, row.names = 1, stringsAsFactors = FALSE)

# transITS data
t.phylum <- read.delim(file = paste0(FILEDIRECTORY, "transITS/Euk_OTU_phylum.txt"), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
t.genus <- read.delim(file = paste0(FILEDIRECTORY, "transITS/Euk_OTU_genus.txt"), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)



##### ---------- DEFINE PARAMETERS

PATIENTS <- colnames(v.genus)
TAXLEVELS <- c("genus","phylum")
MINREAD <- c(2, seq(5, 100, by = 5), seq(200, 500, by = 100))



##### ---------- CALCULATE HAMMING DISTANCE OVER MULTIPLE THRESHOLDS

### Loop over taxonomic ranks
for(i in 1:length(TAXLEVELS)) {

  # Load OTU data
  DB.v <- get(x = paste0("v.", TAXLEVELS[i]), envir = globalenv())
  DB.t <- get(x = paste0("t.", TAXLEVELS[i]), envir = globalenv())

  DB.v <- DB.v[,PATIENTS]
  DB.t <- DB.t[,PATIENTS]

  colnames(DB.v) <- paste0("v.", colnames(DB.v))
  colnames(DB.t) <- paste0("t.", colnames(DB.t))

  # Create result matrix
  mat <- as.data.frame(matrix(nrow = (length(PATIENTS)+1), ncol=2*length(MINREAD)))
  rownames(mat) <- c(PATIENTS,"MaxPossibleScore")
  colnames(mat) <- c(MINREAD, paste0("Bg.", MINREAD)) # "Bg." columns = random background scores
  
  
  ### Loop over minimum read thresholds
  for(j in 1:length(MINREAD)) {

    # Prefilter combined 18S/transITS OTU table by min reads, and convert to binary
    df<- merge(DB.v, DB.t, by = "row.names", all = TRUE)
    df[is.na(df)] <- 0
    rownames(df) <- df$Row.names
    df$Row.names <- NULL
    df[df < MINREAD[j]] <- 0
    df[df >= MINREAD[j]] <- 1
    df <- df[rowSums(df)>0,]

    
    ### Calculate Hamming distance
    df.hamm <- as.data.frame(hamming.distance(as.matrix(t(df))))
    
    # Subset only scores between 18S and transITS
    df.hamm <- df.hamm[grep("v.", rownames(df.hamm)),grep("t.", rownames(df.hamm))]
    rownames(df.hamm) <- gsub("v.", "", rownames(df.hamm))
    colnames(df.hamm) <- gsub("t.", "", colnames(df.hamm))

    # Extract "within-patient" scores
    for(k in 1:length(PATIENTS)) {
      mat[PATIENTS[k],j] <- df.hamm[PATIENTS[k],PATIENTS[k]]
      mat[nrow(mat),j] <- nrow(df)
    }
    
    # Extract "between-patient" scores, averaged over all other patient combinations
    if(all(rownames(df.hamm)==colnames(df.hamm))&all(rownames(df.hamm)==rownames(mat)[1:nrow(mat)-1])) {
      for(k in 1:length(PATIENTS)) {
        mat[k,length(MINREAD)+j] <- mean(as.numeric(df.hamm[k,-k])) 
        mat[nrow(mat),length(MINREAD)+j] <- nrow(df)
      }
    }
  }

  
  ### Save score matrix to global environment
  assign(paste0(TAXLEVELS[i],".hammscore"), data.frame(mat))

  
  ### Write hamming distance output for taxonomic rank
  write.csv(mat, file = paste0(OUTDIR, "HammingSamples_18S.transITS_", TAXLEVELS[i], "_",Sys.Date(), ".csv"), row.names = TRUE)
}

rm(df,df.hamm,DB.t,DB.v,mat)



##### ---------- RUN STATISTICAL TEST AND PLOT HISTOGRAMS


### Define parameters

BIN <- c(0.02,0.03,0.04)

matMEANS.stat <- as.data.frame(matrix(nrow=4, ncol=length(TAXLEVELS))) 
rownames(matMEANS.stat) <- c("method","D.val","p.val","alternativehypothesis")
colnames(matMEANS.stat) <- TAXLEVELS


### Normalize scores, run statistics and plots

for(x in 1:length(TAXLEVELS)) {
  mat <- get(paste0(TAXLEVELS[x],".hammscore"), globalenv())

  # Normalize Hamming distance by maximum scores at each minimum read (MINREAD) cutoff
  # Note: the number of taxa changes with MINREAD filtering
  TOT <- mat[nrow(mat),]
  mat <- mat / TOT[col(mat)]

  # Average normalized scores and write matrix
  matMEANS <- data.frame(means=rowMeans(mat[1:(nrow(mat)-1),1:length(MINREAD)]), meansBg=rowMeans(mat[1:(nrow(mat)-1),(length(MINREAD)+1):ncol(mat)]))
  write.csv(matMEANS, file = paste0(OUTDIR, "HammingSamples_18S.transITS_AveragedScores_", TAXLEVELS[x], "_",Sys.Date(), ".csv"), row.names = TRUE)
  
  MAXSCORE <- ceiling(max(matMEANS))
  
  
  ### Statistical evaluation - Kolmogorov-Smirnov test to compare distributions (nonparametric)

  ks.stat <- ks.test(matMEANS$means, matMEANS$meansBg, alternative = "two.sided")
  matMEANS.stat["method",x] <- ks.stat$method
  matMEANS.stat["D.val",x] <- ks.stat$statistic
  matMEANS.stat["p.val",x] <- ks.stat$p.value
  matMEANS.stat["alternativehypothesis",x] <- ks.stat$alternative

  
  ### Plot histograms 
  
  # Plot using different bin sizes
  for(m in 1:length(BIN)) {
    png(file = paste0(OUTDIR, "HammingSamples_", TAXLEVELS[x], "_Avg_MinRead_2-500_Scores_bin",BIN[m],".png"), width = 2000, height = 1500, res = 300)
    hist(as.matrix(matMEANS[1:nrow(matMEANS), 1]), breaks=seq(0,0.6,by=BIN[m]), 
           col=rgb(1,0,0,0.5), lty="blank", main=paste0(TAXLEVELS[x], " avg scores 2-500 min reads"), 
           xlab = "Normalized hamming distance", xlim = c(0,0.6), ylim=c(0,30))
    hist(as.matrix(matMEANS[1:nrow(matMEANS), 2]), breaks=seq(0,0.6,by=BIN[m]), 
           col=rgb(0.4,0.4,0.8,0.5), lty="blank", add=T)
    box()
    dev.off()
  }
}


### Write statistics table

write.csv(matMEANS.stat, file=paste0(OUTDIR, "ks.test_", TAXLEVELS[x], "_Avg_MinRead_2-500.csv"), row.names = TRUE)



##### ---------- SAVE DATA

save.image(paste0(OUTDIR,"Hamming_transITS18S_associations.rdata"))


# A. Popovic, Parkinson Lab. July 2018