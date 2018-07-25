#################################################################
#####       OTU associations / ordination using sPLS-DA      #####



##### ---------- LOAD LIBRARIES, SET DIRECTORIES AND IMPORT DATA

library("mixOmics")
library("metagenomeSeq")

dir.create("OTU_sPLSDA")
WD <- paste0(getwd(), "/", "OTU_sPLSDA/")
setwd(WD)

FILEPATH <- "C:Users/User/ProjectData/"


### Import OTU tables

# columns are samples, rows are OTUs
eukOTU <- read.delim(paste0(FILEPATH, "Euk_genusOTU.txt"), sep = "\t", stringsAsFactors = FALSE)
bacOTU <- read.delim(paste0(FILEPATH, "Bac_genusOTU.txt"), sep = "\t", stringsAsFactors = FALSE)

OTU <- bacOTU


### Import OTU taxonomies and sample metadata

# OTU taxonomy: columns are taxonomic ranks, rows are OTUs
# sample metadata: columns are various measures/features of interest, rows are samples
taxa_bac <- read.delim(paste0(FILEPATH, "OTU_taxonomies.txt"), header=TRUE, row.names = 1, stringsAsFactors = FALSE)
sampledata <- read.delim(paste0(FILEPATH, "sample_metadata.txt"), sep = "\t", row.names = 1, stringsAsFactors = FALSE)

sampledata$TotalEukReads <- colSums(eukOTU[,rownames(sampledata)])


### Select/prepare variable(s) for ordination

# Blastocystis, presence >5 reads
sampledata$Blastocystis <- t(eukOTU["Blastocystis", rownames(sampledata)])
sampledata$Blastocystis.r <- sampledata$Blastocystis / sampledata$TotalEukReads

sampledata$Blast.min5reads <- sampledata$Blastocystis
sampledata[sampledata$Blastocystis<5,"Blast.min5reads"] <- "N"
sampledata[sampledata$Blastocystis>=5,"Blast.min5reads"] <- "Y"

variables <- sampledata$Blast.min5reads



##### ---------- FILTER OTUs

### Filter OTUs 

# >=5 reads in at least 10% of samples
OTU <- OTU[,rownames(variables)] #reorder columns
OTU.filtered <- OTU[rowSums(OTU>=5)>=ceiling(ncol(OTU)*0.1),]



##### ---------- NORMALIZE OTU DATA

### Cumulative Sum Scaling

# Create MetagenomeSeq object
OBJ.variables <- AnnotatedDataFrame(data = as.data.frame(variables)) # Sample data
OBJ <- newMRexperiment(OTU.filtered, phenoData=OBJ.variables)

# Normalize OTU data
p <- cumNormStat(OBJ)
OBJ.CSS <- cumNorm(OBJ, p=p)
OTU.CSS <- t(MRcounts(OBJ.CSS, norm=TRUE, log=TRUE)) #save data 
rm(OBJ.variables,OBJ,p,OBJ.CSS)



##### ---------- Preliminary PCA with 10 components

df <- OTU.CSS

df.pca <- pca(df, ncomp = 10, center = TRUE, scale = TRUE)

plotIndiv(df.pca, comp = c(1,2), group = as.factor(variables), ind.names = FALSE, ellipse = TRUE, ellipse.level = 0.95, legend = TRUE, 
          centroid = TRUE, style="ggplot2", title = 'PCA Comp 1-2', size.xlabel = rel(1.2), 
          size.ylabel = rel(1.2), size.legend.title = rel(1.1), legend.title = "Blastocystis",
          size.axis = rel(1.2), pch = c(19,17), col = c('blue','red'), legend.position = 'right')

plot(df.pca)
# Not great separation in PCA. (s)PLS-DA anlaysis to explore clusters.



##### ---------- (s)PLS-DA ON OTU DATA

### Initial calculation and plot

df.plsda <- plsda(df, as.factor(variables), ncomp = 10)

pls.plot <- plotIndiv(df.plsda, ind.names = FALSE, ellipse = TRUE, ellipse.level = 0.95, legend = TRUE, 
          centroid = TRUE, style="ggplot2", title = 'PLS-DA Comp 1-2', size.xlabel = rel(1.2), 
          size.ylabel = rel(1.2), size.legend.title = rel(1.1), legend.title = "Blastocystis",
          size.axis = rel(1.2), pch = c(19,17), col = c('lightblue','red2'), legend.position = 'right')

plotLoadings(df.plsda, comp = 1, method = 'mean', contrib = 'max', ndisplay = 30)


### Evaluate initial model

# Use 5-fold cross-validation, repeated 50x, to identify optimal number of components
#set.seed(2543)
perf.plsda <- perf(df.plsda, validation = "Mfold", folds = 10, 
                   progressBar = TRUE, auc = TRUE, nrepeat = 50) 

perf.plsda$error.rate  # error rates

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")


### sPLS-DA - tune to select number of components and features

# create list of values tested for each component
#list.keepX <- c(1:4,seq(5, 50, 5)) 
list.keepX <- c(seq(5, 50, 5)) 

#set.seed(2543)
tune.splsda <- tune.splsda(df, as.factor(variables), ncomp = 10, validation = 'Mfold', folds = 10, 
                           progressBar = TRUE, dist = 'max.dist',
                           test.keepX = list.keepX, nrepeat = 50)

# view optimal number of components, and features per component
tune.splsda$error.rate
tune.splsda$choice.ncomp$ncomp
tune.splsda$choice.keepX
plot(tune.splsda, col = color.jet(10))

choice.ncomp <- 2 #tune.splsda$choice.ncomp$ncomp
choice.keepX <- tune.splsda$choice.keepX[1:choice.ncomp] 


### sPLS-DA with tuned parameters

df.splsda <- splsda(df, as.factor(variables), ncomp = choice.ncomp, keepX = choice.keepX)

pls.plot <- plotIndiv(df.splsda, ind.names = FALSE, ellipse = TRUE, ellipse.level = 0.95, legend = TRUE, 
            centroid = TRUE, style="ggplot2", title = 'PLS-DA Comp 1-2', size.xlabel = rel(1.2), 
            size.ylabel = rel(1.2), size.legend.title = rel(1.1), legend.title = "Blastocystis",
            size.axis = rel(1.2), pch = c(19,17), col = c('lightblue','red2'), legend.position = 'right')
            #Ellipses: confidence level set to 0.95 per default

plotLoadings(df.splsda, comp = 1, method = 'mean', contrib = 'max')


### Evaluate errors through 10-fold cross-validation

perf.splsda <- perf(df.splsda, validation = "Mfold", folds = 10, 
                    progressBar = TRUE, auc = TRUE, nrepeat = 50) 

plot(perf.splsda, sd = TRUE, legend.position = "horizontal")

perf.splsda$error.rate
perf.splsda$error.rate.class


# Output final selection of features with their weight coefficient
selectVar(df.splsda, comp = 1)$value
plot(perf.splsda$features$stable[[1]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 1', las =2)

selectVar(df.splsda, comp = 2)$value
plot(perf.splsda$features$stable[[2]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 2', las =2)


# Leave one out error calculation
perf.splsda.loo <- perf(df.splsda, validation = "loo", criterion = "all", progressBar = TRUE) 
perf.splsda.loo$error.rate.all


### Obtain R2 and Q2 values for sPLS-DA

# (Not produced for sPLS-DA by perf function)
X <- df
Y <- as.factor(variables)
Y.mat <- unmap(Y)
res <- spls(X,Y.mat, ncomp = choice.ncomp, keepX = choice.keepX)
val <- perf(res, criterion = c("R2", "Q2"), validation = "Mfold", folds = 10, 
            progressBar = TRUE, nrepeat = 50)
val
val$Q2.total
val$R2


### Proportion of variance explained by PLS components
Rd.YvsU <- cor(as.numeric(as.factor(variables)), df.splsda$variates$X[, 1:choice.ncomp])
Rd.YvsU <- apply(Rd.YvsU^2, 2, sum)
Rd.Y <- cbind(Rd.YvsU, cumsum(Rd.YvsU))
colnames(Rd.Y) <- c("Proportion", "Cumulative")

# View proportion and cumulative values
Rd.Y



##### ---------- SAVE DATA

save.image(paste0(WD,"sPLSDA.rdata"))


# A. Popovic, Parkinson Lab. July 2018