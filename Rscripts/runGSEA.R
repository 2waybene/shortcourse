
gct.file <- "U:/My Documents/NIEHS_shortCourse/summer2016/GSEA/DEG_overall_gender_diff.txt"
df <- read.delim(gct.file , header=TRUE, sep="\t", skip=2,  blank.lines.skip=TRUE)

source("U:/My Documents/NIEHS_shortCourse/Rscripts/GSEA_1_0.R")
setwd("U:/My Documents/NIEHS_shortCourse/summer2016/")

GSEA(
  # Input/Output Files :-------------------------------------------
  input.ds = "GSEA/dt_no_duplicate.gct",
  input.cls = "GSEA/gender_label_word.cls",
  gs.db = "GSEA/c4.cm.v5.1.symbols.gmt",
  output.directory = "GSEA_OUT/",
  # Program parameters :-----------------------------------------
  doc.string = "Gender_C1",
  non.interactive.run = FALSE,
  reshuffling.type = "sample.labels",
  nperm = 1000,
  weighted.score.type = 1,
  nom.p.val.threshold = -1,
  fwer.p.val.threshold = -1,
  fdr.q.val.threshold = 0.25,
  topgs = 20,
  adjust.FDR.q.val = FALSE,
  gs.size.threshold.min = 15,
  gs.size.threshold.max = 500,
  reverse.sign = FALSE,
  preproc.type = 0, 
  random.seed = 111,
  # Tweaks for experts only :-----------------------------------------
  perm.type = 0,
  fraction = 1.0,
  replace = FALSE,
  save.intermediate.results = FALSE,
  OLD.GSEA = FALSE,
  use.fast.enrichment.routine = TRUE
)


