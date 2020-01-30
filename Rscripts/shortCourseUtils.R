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



affyProb2Symbols.hu <- function (aListOfProbes, db = hgu133plus2.db)
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

