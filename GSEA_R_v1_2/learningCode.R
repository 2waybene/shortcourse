
##===========================================================
##  This is different version of implementation
##===========================================================
library(RGSEA)
if(interactive()) {
  data(e1)
  data(e2)
  RGSEAfix(e1,e2, queryclasses=colnames(e1), refclasses=colnames(e2),
           random=20000, featurenum=1000, iteration=100)->test
}

str(e1)
str(e2)

##===========================================================
##  This is GSEA 1.2 from BroadInstitute
##===========================================================


library("utils")
library("tools")
library("dplyr")
library("GSEA")

##  simple example and it works

GSEA(input.ds = system.file('extdata', 'Leukemia_hgu95av2.gct', package = 'GSEA', mustWork = TRUE),
     input.cls = system.file('extdata', 'Leukemia.cls', package = 'GSEA', mustWork = TRUE),
     input.chip = system.file('extdata', 'Human_AFFY_HG_U95_MSigDB_7_0_final.chip',
                              package = 'GSEA', mustWork = TRUE), gs.db = system.file('extdata',
                                                                                      'h.all.v7.0.symbols.gmt', package = 'GSEA', mustWork = TRUE),
     collapse.dataset = TRUE, collapse.mode = 'max')




##===============================
##  does not work here
##===============================

input.ds = "C:/Users/li11/Documents/R/win-library/3.6/GSEA/extdata/Leukemia_hgu95av2.gct"                     # Input gene expression dataset file in GCT format
input.cls = "C:/Users/li11/Documents/R/win-library/3.6/GSEA/extdata/Leukemia.cls"                             # Input class vector (phenotype) file in CLS format
gs.db = "C:/Users/li11/Documents/R/win-library/3.6/GSEA/extdata/h.all.v7.0.symbols.gmt"                       # Gene set database in GMT format
input.chip = "C:/Users/li11/Documents/R/win-library/3.6/GSEA/extdata/Human_AFFY_HG_U95_MSigDB_7_0_final.chip" # CHIP File
output.directory      = "X:/myGit/shortcourse/GSEA_R_v1_2/outputDir/"        # Directory where to store output and results (default: "")

GSEA(
  # Input/Output Files :-------------------------------------------------------------------------------
  #input.ds = inputds,                    # Input gene expression dataset file in GCT format
  #input.cls = inputcls,                  # Input class vector (phenotype) file in CLS format
  #gs.db = gsdb,                          # Gene set database in GMT format
  #input.chip = inputchip,               # CHIP File
  #output.directory      = outdir,        # Directory where to store output and results (default: "")
  
  input.ds = "C:/Users/li11/Documents/R/win-library/3.6/GSEA/extdata/Leukemia_hgu95av2.gct"    ,                 # Input gene expression dataset file in GCT format
  input.cls = "C:/Users/li11/Documents/R/win-library/3.6/GSEA/extdata/Leukemia.cls"             ,                # Input class vector (phenotype) file in CLS format
  gs.db = "C:/Users/li11/Documents/R/win-library/3.6/GSEA/extdata/h.all.v7.0.symbols.gmt"        ,               # Gene set database in GMT format
  input.chip = "C:/Users/li11/Documents/R/win-library/3.6/GSEA/extdata/Human_AFFY_HG_U95_MSigDB_7_0_final.chip", # CHIP File
  output.directory      = "X:/myGit/shortcourse/GSEA_R_v1_2/outputDir/"     ,   # Directory where to store output and results (default: "")
  
  #  Program parameters :-------------------------------------------------------------------------------
#  doc.string            = outname,         # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
#  non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
#  reshuffling.type      = permutation,     # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels"
#  nperm                 = as.integer(nperms),            # Number of random permutations (default: 1000)
  weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
  nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
  fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
  fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
  topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
  adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
 # gs.size.threshold.min = as.integer(minsize),         # Minimum size (in genes) for database gene sets to be considered (default: 25)
#  gs.size.threshold.max = as.integer(maxsize),         # Maximum size (in genes) for database gene sets to be considered (default: 500)
  reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
  preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
 # random.seed           = as.integer(as.POSIXct(Sys.time())),            # Random number generator seed. (default: 123456)
  perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
  fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
  replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
  collapse.dataset      = T, #collapsedataset, # Collapse dataset to gene symbols using a user provided chip file (default: F)
  collapse.mode         = 2, #collapsemode,
  save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
  use.fast.enrichment.routine = T,         # Use faster routine to compute enrichment for random permutations (default: T)
#  gsea.type = rankmethod,                     # Select Standard GSEA (default) or preranked
#  rank.metric = rankmetric
)

GSEA(input.ds = system.file('extdata', 'Leukemia_hgu95av2.gct', package = 'GSEA', mustWork = TRUE),
     input.cls = system.file('extdata', 'Leukemia.cls', package = 'GSEA', mustWork = TRUE),
     input.chip = system.file('extdata', 'Human_AFFY_HG_U95_MSigDB_7_0_final.chip',
                              package = 'GSEA', mustWork = TRUE), gs.db = system.file('extdata',
                                                                                      'h.all.v7.0.symbols.gmt', package = 'GSEA', mustWork = TRUE),
     collapse.dataset = TRUE, collapse.mode = 'max')

