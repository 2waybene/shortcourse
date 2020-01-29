


library(affy);
library(limma);
library(affycoretools)
library(annotate)
library(mouse4302.db)


##=================================================
##  functions
##=================================================

affyProb2Symbols <- function (aListOfProbes, db = mouse4302.db)
{
  require(magrittr)
  require(mouse4302.db)
  probeAnnoted <- AnnotationDbi::select(
    x       = db,
    keys    = aListOfProbes,
    columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
    keytype = "PROBEID"
  )
  return (probeAnnoted )
}


##=========================================
##  Learning scratch starts
##=========================================
library(GEOquery)
GEOset <- getGEO("GSE40661") # returns a list of 'ExpressionSet's
GEOset <- GEOset[[1]] 
dim(GEOset)
featureData(GEOset)



##=============================================================================
##  Download mouse430 v2 probes
##  https://www.ebi.ac.uk/arrayexpress/files/A-GEOD-17109/A-GEOD-17109.adf.txt
##=============================================================================


##  testing probes to symbole with the following modules

##  probes directly from database
test.id <- head(keys(mouse4302.db),20)

##  probes directly download from ebi
id.in <- read.table ("x:/project2020/SEM/KEGG_analysis/A-GEOD-17109.adf.txt", header = TRUE, sep="\t")
str(id.in)
probes <- as.character(id.in$AffyProbe)
test.id <- probes[1:20]


##  probes directly from KEGG analysis
kegg.id <- tolower(probe.ids[[1]])

gsub(" ", "", kegg.id) %in% probes
gsub(" ", "", kegg.id) %in% keys(mouse4302.db)

affyProb2Symbols(gsub(" ", "", kegg.id))

##=========================================
##  Learning scratch ends
##=========================================


##==============================================

kegg.results <- read.csv("x:/project2019/SEM/IPA/KEGG_analysis/kegg_28_pathways.csv", header = TRUE)

dim(kegg.results)
str(kegg.results)
head(kegg.results)

probe.ids <- strsplit (as.character(kegg.results$Genes[1]), ",")

require(doParallel)
numCores <- detectCores()/2
cl <- makeCluster(numCores , type='PSOCK')
registerDoParallel(cl)

rows <- dim (kegg.results)[1]
ptime <- system.time({
  r <- foreach(i = 1:rows, .combine='c') %dopar% {
  ##  r <- foreach(i = 1:rows, .combine='cbind') %dopar% {  ## does NOT work here!!
    id <- as.character(kegg.results$KEGGIndex[i])
    probe.ids <- strsplit (as.character(kegg.results$Genes[i]), ",")
    kegg.id <- tolower(probe.ids[[1]])
    geneSymbols <- affyProb2Symbols(gsub(" ", "", kegg.id))$SYMBOL
    id.coverted <- as.data.frame(geneSymbols)
    colnames(id.coverted) <- id
    id.coverted
  }
})[3]
print (paste ("This is parallel time: " , ptime, sep=""))
stopCluster(cl)  ## Modify by JYL, 10/23/19

##===================================
##  why cbind does NOT work??
##  why geneSymbols is factor??
##===================================

##  write out coverted gene symbols for 
##  bootstrap analysis in Rshiny app

dim(kegg.results)
keggSymbols <- c()
for (i in 1:length(r))
{
  keggSymbols [i] <- toupper(paste (r[[i]], collapse = ","))
}
kegg.out <- cbind (kegg.results[,-5], as.data.frame(keggSymbols))

write.table(kegg.out,file = "x:/project2020/SEM/KEGG_analysis/KEGG_4_bootstrap.txt", sep = "\t", row.name = FALSE)

library("xlsx")
write.xlsx(kegg.out, file = "x:/project2020/SEM/KEGG_analysis/test.xlsx", row.names = FALSE, append = FALSE)
##=======================================================


##============================================================================
##  Finished bootstrap analysis with two folders containing the results
##  x:/project2020/SEM/KEGG_analysis/withoutReplacement
##  x:/project2020/SEM/KEGG_analysis/withReplacement
##  Now, getting the analysis for the index-like test
##  On linux desktop: 01934770:
##  calculate_index_streamline_linux_w_arguments.R
##=============================================================================


##  Usage: Rscript --vanilla calculate_index_streamline_linux_w_arguments.R tScoreData path2WithReplacement  path2WithoutReplacement helperScripts  outputDir

##   /usr/lib64/R/bin/exec/R --slave --no-restore --vanilla --file=/ddn/gs1/home/li11/ddnDrive/project2019/SEM/rscripts/calculate_index_streamline_linux_w_arguments.R 
##  --args /ddn/gs1/home/li11/ddnDrive/project2019/SEM/data/new_t_scores_w_lev_03272019.rda /ddn/gs1/  /home/li11/ddnDrive/project2020/SEM/KEGG_analysis/withReplacement/ 
##  /ddn/gs1/home/li11/ddnDrive/project2020/SEM/KEGG_analysis/withoutReplacement/ /ddn/gs1/home/li11/ddnDrive/project2019/SEM/rscripts/helpScripts/SEM_util.R 
##  /ddn/gs1/home/li11/ddnDrive/project2020/SEM/KEGG_analysis/bootstrapIndexResults/
  





