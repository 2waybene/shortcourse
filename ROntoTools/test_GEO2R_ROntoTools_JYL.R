
##==============================================================
##  File Name : geo2R_lungCancer_GSA.R
##  Author    : Jianying Li
##  Comment   : Boxplot for selected GEO samples "GSE18842"
##==============================================================
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE18842", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("10101101010101011010101010101010010110101010101101",
               "01000001010101010101010010110101110101010")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

setwd("X:/myGit/shortcourse/workWgeo/")
write.table(tT, file="geo2r_limma_example.txt", row.names=F, sep="\t")

head(tT)
dim(tT)


limma.all <- topTable(fit2, adjust="fdr", sort.by="B", number=40000)
limma.all <- subset(limma.all, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
dim(limma.all)
length(which(limma.all$adj.P.Val <= 0.05))
GSE18842.limma.all = limma.all[which(limma.all$adj.P.Val <= 0.05),]
str(GSE18842.limma.all)

geneAnnoated <- affyProb2Symbols(tolower(GSE18842.limma.all$ID))
library(tidyverse)
filtered = geneAnnoated %>% distinct(PROBEID, .keep_all = TRUE)
head(filtered)
dim(filtered)


##=====================================
##  Prepare files for GSA now
##=====================================

##  Annotate the affy probes

affyProb2Symbols <- function (aListOfProbes, db = hgu133plus2.db)
{
  require(magrittr)
  require(hgu133plus2.db)
  probeAnnoted <- AnnotationDbi::select(
    x       = db,
    keys    = aListOfProbes,
    columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
    keytype = "PROBEID"
  )
  return (probeAnnoted )
}



geneAnnoated <- affyProb2Symbols(tolower(rownames(ex)))
library(tidyverse)
filtered = geneAnnoated %>% distinct(PROBEID, .keep_all = TRUE)
head(filtered)


data <- list(x = as.matrix(ex), y = y, geneid = paste ("g", 1:dim(ex)[1], sep=""), genenames = filtered$SYMBOL, logged2 = TRUE)


##=====================================
##  Run limma analysis on this
##=====================================
library(limma)

# counts is the assay data I already have.
dim(ex)
# [1] 64102     8

# Creates a new ExpressionSet object (quite bare, only the assay data)
asdf <- ExpressionSet(assayData = ex)

# Returns the data you put in.
exprs(asdf)



##=======================================================================
load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))
head(top)
dim(top)

entrez = paste ("hsa:", filtered$ENTREZID, sep="")
dim(GSE18842.limma.all)
dat <- cbind(GSE18842.limma.all[,c(6,3,2)], entrez)
dat$entrez <- as.character(entrez)
dat <- dat [ -c(which(dat$entrez == "hsa:NA")),]
top <- dat
dim(dat)
