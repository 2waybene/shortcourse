##===========================================================
##  NSCLC_used_in_Nguyen_2019.R
##  GEO_ID	Disease	Normal	Condition	Pubmed_ID	Tissue	Platform
##  GSE18842	Non-small cell lung cancer	44	44	20878980	Lung      	HG-U133 Plus 2.0
##  GSE19188	Non-small cell lung cancer	62	91	20421987	Lung      	HG-U133 Plus 2.0
##  GSE19804	Non-small cell lung cancer	60	60	20802022	Lung      	HG-U133 Plus 2.0
##  GSE50627	Non-small cell lung cancer	6	9	25881239	Lung      	HuGene-10st
##  GSE6044	Non-small cell lung cancer	5	31	18992152	Lung      	Hu-Focus

##===========================================================
##  GSE18842
##  Differential expression analysis with limma
##===========================================================
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

#tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=dim(ex)[1])
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file="X:/myGit/shortcourse/ROntoTools/NSCLC_data/GSE18842_lmFit_all.txt", row.names=F, sep="\t")
save(tT, file="X:/myGit/shortcourse/ROntoTools/NSCLC_data/GSE18842_lmFit_all.rda")


##===========================================================
##  GSE19188
##  Differential expression analysis with limma
##===========================================================
################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE19188", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("10000110100010101010110101010101010110101010110100",
               "10111011010101101010011011010101010100100110010100",
               "10111111010011101110111110101101110111111101011111",
               "101010")
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
#tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)


tT <- topTable(fit2, adjust="fdr", sort.by="B", number=dim(ex)[1])
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file="X:/myGit/shortcourse/ROntoTools/NSCLC_data/GSE19188_lmFit_all.txt", row.names=F, sep="\t")
save(tT, file="X:/myGit/shortcourse/ROntoTools/NSCLC_data/GSE19188_lmFit_all.rda")


##===========================================================
##  GSE19804
##  Differential expression analysis with limma
##===========================================================
################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE19804", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "undefined"
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
cont.matrix <- makeContrasts(Gu-Gd, Ge-Gd, Gf-Ge, Gi-Gf, Gn-Gi, Gu-Gn, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
#tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)


tT <- topTable(fit2, adjust="fdr", sort.by="B", number=dim(ex)[1])
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file="X:/myGit/shortcourse/ROntoTools/NSCLC_data/GSE19804_lmFit_all.txt", row.names=F, sep="\t")
save(tT, file="X:/myGit/shortcourse/ROntoTools/NSCLC_data/GSE19804_lmFit_all.rda")


##===========================================================
##  GSE50627
##  Differential expression analysis with limma
##===========================================================

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE50627", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "undefined"
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
cont.matrix <- makeContrasts(Gu-Gd, Ge-Gd, Gf-Ge, Gi-Gf, Gn-Gi, Gu-Gn, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
#tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

#tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
#write.table(tT, file=stdout(), row.names=F, sep="\t")

tT <- topTable(fit2, adjust="fdr", sort.by="B", number=dim(ex)[1])
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file="X:/myGit/shortcourse/ROntoTools/NSCLC_data/GSE50627_lmFit_all.txt", row.names=F, sep="\t")
save(tT, file="X:/myGit/shortcourse/ROntoTools/NSCLC_data/GSE50627_lmFit_all.rda")


##===========================================================
##  GSE6044
##  Differential expression analysis with limma
##===========================================================
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE6044", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL201", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "undefined"
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
cont.matrix <- makeContrasts(Gu-Gd, Ge-Gd, Gf-Ge, Gi-Gf, Gn-Gi, Gu-Gn, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
#tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

#tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
#write.table(tT, file=stdout(), row.names=F, sep="\t")


tT <- topTable(fit2, adjust="fdr", sort.by="B", number=dim(ex)[1])
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file="X:/myGit/shortcourse/ROntoTools/NSCLC_data/GSE50627_lmFit_all.txt", row.names=F, sep="\t")
save(tT, file="X:/myGit/shortcourse/ROntoTools/NSCLC_data/GSE50627_lmFit_all.rda")


