
##===================================================
##  File Name : GEO2R_lungCancer_ROntoTools.R
##  Author    : Jianying Li
##===================================================


##=============================================================
##  Part 1 :  Download a selected GEO samples "GSE18842"
##            apply limma.fit
##==============================================================


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


geneAnnoated <- affyProb2Symbols(tolower(GSE18842.limma.all$ID))
library(tidyverse)
filtered = geneAnnoated %>% distinct(PROBEID, .keep_all = TRUE)

#head(filtered)
entrez = paste ("hsa:", filtered$ENTREZID, sep="")
#dim(GSE18842.limma.all)
dat <- cbind(GSE18842.limma.all[,c(6,3,2)], entrez)
dat$entrez <- as.character(entrez)
dat <- dat [ -c(which(dat$entrez == "hsa:NA")),]
top <- dat
#dim(dat)

# data <- list(x = as.matrix(ex), y = y, geneid = paste ("g", 1:dim(ex)[1], sep=""), genenames = filtered$SYMBOL, logged2 = TRUE)


##=======================================================================







##=============================================================
##  Part 2 :  Running ROntoTools
##            Code from rontotools
##==============================================================

###################################################
### code chunk number 1: rontotools.Rnw:55-58
###################################################
require(graph)
require(ROntoTools)
kpg <- keggPathwayGraphs("hsa", verbose = FALSE)


###################################################
### code chunk number 2: rontotools.Rnw:62-63 (eval = FALSE)
###################################################
## kpg <- keggPathwayGraphs("hsa", updateCache = TRUE, verbose = TRUE)

#length(keggPathwayNames(organism = "hsa", updateCache = FALSE, verbose = TRUE))
#[1] 337

###################################################
### code chunk number 3: rontotools.Rnw:69-70
###################################################
#head(names(kpg))
#length(kpg)
#[1] 221

###################################################
### code chunk number 4: rontotools.Rnw:74-77
###################################################
#kpg[["path:hsa04110"]]
#head(nodes(kpg[["path:hsa04110"]]))
#head(edges(kpg[["path:hsa04110"]]))


###################################################
### code chunk number 5: rontotools.Rnw:81-82
###################################################
#head(edgeData(kpg[["path:hsa04110"]], attr = "subtype"))


###################################################
### code chunk number 6: rontotools.Rnw:86-90
###################################################
kpg <- setEdgeWeights(kpg, edgeTypeAttr = "subtype", 
                      edgeWeightByType = list(activation = 1, inhibition = -1, 
                                              expression = 1, repression = -1), 
                      defaultWeight = 0)


###################################################
### code chunk number 7: rontotools.Rnw:94-95
###################################################
#head(edgeData(kpg[["path:hsa04110"]], attr = "weight"))


###################################################
### code chunk number 8: rontotools.Rnw:99-101
###################################################
kpn <- keggPathwayNames("hsa")
#head(kpn)

#length(kpn)
#[1] 337
#kpn[which(names(kpn) %in% "path:hsa05216")]

###################################################
### code chunk number 10: rontotools.Rnw:127-136
###################################################
fc <- top$logFC[top$adj.P.Val <= .01]
names(fc) <- top$entrez[top$adj.P.Val <= .01]

pv <- top$P.Value[top$adj.P.Val <= .01]
names(pv) <- top$entrez[top$adj.P.Val <= .01]

#head(fc)
#head(pv)


###################################################
### code chunk number 11: rontotools.Rnw:140-145
###################################################
fcAll <- top$logFC
names(fcAll) <- top$entrez

pvAll <- top$P.Value
names(pvAll) <- top$entrez


###################################################
### code chunk number 12: rontotools.Rnw:149-151
###################################################
ref <- top$entrez
#head(ref)


###################################################
### code chunk number 13: rontotools.Rnw:163-165
###################################################
kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)
#head(nodeWeights(kpg[["path:hsa04110"]]))


###################################################
### code chunk number 14: rontotools.Rnw:177-178
###################################################
peRes <- pe(x = fc, graphs = kpg, ref = ref,  nboot = 200, verbose = FALSE)


###################################################
### code chunk number 15: rontotools.Rnw:182-186
###################################################
head(Summary(peRes))

head(Summary(peRes, pathNames = kpn, totalAcc = FALSE, totalPert = FALSE,
             pAcc = FALSE, pORA = FALSE, comb.pv = NULL, order.by = "pPert"))

peSummary <- Summary(peRes, pathNames = kpn, totalAcc = FALSE, totalPert = FALSE,
                     pAcc = FALSE, pORA = FALSE, comb.pv = NULL, order.by = "pPert")

write.table (list (KEGG = peSummary$pathNames), file = "x:/project2019/RNAseqProj/ESCC/lungCancer/pathway.txt", sep = "\t", row.names = FALSE)
peSummary[13,]


###################################################
### code chunk number 16: peRes_twoway1
###################################################
plot(peRes)


###################################################
### code chunk number 17: peRes_twoway2
###################################################
plot(peRes, c("pAcc", "pORA"), comb.pv.func = compute.normalInv, threshold = .01)


###################################################
### code chunk number 20: pePathway_twoway_Acc
###################################################
plot(peRes@pathways[["path:hsa05216"]], type = "two.way")


###################################################
### code chunk number 21: pePathway_boot_Acc
###################################################
plot(peRes@pathways[["path:hsa05216"]], type = "boot")


###################################################
### code chunk number 22: fig3
###################################################
plot(peRes@pathways[["path:hsa05216"]], type = "two.way")



###################################################
### code chunk number 24: pePathway_graph_Pert
###################################################
p <- peRes@pathways[["path:hsa05216"]]
g <- layoutGraph(p@map, layoutType = "dot")
graphRenderInfo(g) <- list(fixedsize = FALSE)
edgeRenderInfo(g) <- peEdgeRenderInfo(p)
nodeRenderInfo(g) <- peNodeRenderInfo(p)
renderGraph(g)


###################################################
### code chunk number 25: pePathway_graph_Pert2
###################################################
p <- peRes@pathways[["path:hsa04660"]]
g <- layoutGraph(p@map, layoutType = "dot")
graphRenderInfo(g) <- list(fixedsize = FALSE)
edgeRenderInfo(g) <- peEdgeRenderInfo(p)
nodeRenderInfo(g) <- peNodeRenderInfo(p)
renderGraph(g)


###################################################
### p53 signaling pathway
###################################################
p <- peRes@pathways[["path:hsa04115"]]
g <- layoutGraph(p@map, layoutType = "dot")
graphRenderInfo(g) <- list(fixedsize = FALSE)
edgeRenderInfo(g) <- peEdgeRenderInfo(p)
nodeRenderInfo(g) <- peNodeRenderInfo(p)
renderGraph(g)



peSummary[66,]
#pathNames      pPert pPert.fdr
#path:hsa05222 Small cell lung cancer 0.02487562 0.0793825

peSummary[106,]
#pathNames     pPert pPert.fdr
#path:hsa05223 Non-small cell lung cancer 0.1144279 0.2342533
###################################################
### Small cell lung cancer
###################################################
p <- peRes@pathways[["path:hsa05222"]]
g <- layoutGraph(p@map, layoutType = "dot")
graphRenderInfo(g) <- list(fixedsize = FALSE)
edgeRenderInfo(g) <- peEdgeRenderInfo(p)
nodeRenderInfo(g) <- peNodeRenderInfo(p)
renderGraph(g)

###################################################
### Non-small cell lung cancer
###################################################
p <- peRes@pathways[["path:hsa05223"]]
g <- layoutGraph(p@map, layoutType = "dot")
graphRenderInfo(g) <- list(fixedsize = FALSE)
edgeRenderInfo(g) <- peEdgeRenderInfo(p)
nodeRenderInfo(g) <- peNodeRenderInfo(p)
renderGraph(g)

