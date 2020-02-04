##==============================================================
##  File Name : geo2R_lungCancer_GSA.R
##  Author    : Jianying Li
##  Comment   : Boxplot for selected GEO samples "GSE18842"
##==============================================================


##=====================================
##  Download data from GEO
##=====================================
library(Biobase)
library(GEOquery)

# load series and platform data from GEO
gset <- getGEO("GSE18842", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- paste0("10101101010101011010101010101010010110101010101101",
               "01000001010101010101010010110101110101010")

sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  #set group names

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("control","tumor")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE18842", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")


##=====================================
##  Running the SAM analysis
##  Not the main purpose here!
##=====================================

##  getting label

y <- rep ("1", dim(ex)[2])
y[sml %in% "G1"] = 2
#y
samfit <- SAM(as.matrix(ex), y, resp.type = "Two class unpaired")
plot(samfit)


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



##==================================
## Load gene set database
##==================================

clean.gmt.data <- function (db)
{
  # db <- c1 
  # str(db)
  # db$geneset.names
  # db$geneset.descriptions
  temp.geneset <- c(list())
  for (i in 1:length(db$genesets))
  {
    temp.geneset[[i]] <- db$genesets[[i]]
    if (length(which(temp.geneset[[i]]=="")) > 0 )
    {
      temp.geneset[[i]]  <- temp.geneset[[i]][-which(temp.geneset[[i]]=="")]
    }
  }
  db$genesets <- temp.geneset 
  return (db)
}


# original file
c1 <- GSA.read.gmt("X:/myGit/shortcourse/GSA_R/Knowledge_based_gene_sets_20161005.gmt")
#c1 <- GSA.read.gmt("cleaningData/genset_2_correct.txt")
#c1 <- GSA.read.gmt("cleaningData/genset_2_correct_from_1.txt")
#c1 <- GSA.read.gmt("cleaningData/genset_2.txt")
# GSA.read.gmt constructs a list of 3

c1.cleaned <- clean.gmt.data(c1)
#str(c1)
#str(c1$genesets)

print ("Gene sets in database file (20161005).gmt")
c1$geneset.names


##=====================================
##  Running the GSA now
##=====================================

## Setting the random seed for future validation
set.seed(1234567)

c1 <- c1.cleaned

gsafit <- GSA(as.matrix(ex), y, c1$genesets, genenames = filtered$SYMBOL, resp.type = "Two class unpaired", nperms = 1000)

## Do I need the unique gene symbol??
#gsafit <- GSA(x, y, c1$genesets, unique(z$GeneSymbol), resp.type = "Two class unpaired", nperms = 1000)

str(gsafit)
plot(gsafit$GSA.scores)


##=====================================
##  Running the GSA now
##  Second
##=====================================
gsafit$fdr.lo
gsafit$fdr.hi

# To many and will NOT print it to the html file
# gsafit$gene.scores

GSA.plot(gsafit)

GSA.listsets(gsafit, c1$geneset.names, FDRcut=.5)

GSA.correlate(c1$genesets, filtered$SYMBOL)
gsa2 <- GSA.func(as.matrix(ex), y, c1$genesets, filtered$SYMBOL, resp.type = "Two class unpaired")
#gsa2 <- GSA.func(x, y, c1$genesets, unique(z$GeneSymbol), resp.type = "Two class unpaired")
str(gsa2)

#for loop from 1:length(gsals$negative[,1])
GSA.genescores(5, c1$genesets, gsafit,filtered$SYMBOL)
GSA.plot(gsafit, fac = 1, FDRcut = 1)

gsals <- GSA.listsets(gsafit, c1$geneset.names, FDRcut=.5)
str(gsals)
str(gsals$negative)
gsals$negative[,1]
length(gsals$negative[,1])
as.numeric(gsals$negative[,1][1])
GSA.genescores(as.numeric(gsals$negative[,1][1]), c1$genesets, gsafit, filtered$SYMBOL)

