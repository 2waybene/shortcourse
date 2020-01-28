library(RGSEA)
if(interactive()) {
  data(e1)
  data(e2)
  RGSEAfix(e1,e2, queryclasses=colnames(e1), refclasses=colnames(e2),
           random=20000, featurenum=1000, iteration=100)->test
}

str(e1)
str(e2)
