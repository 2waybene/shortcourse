setwd("x:/project2016/woychik-ma-rnaseq/comparison/")
library(rmarkdown)
library(knitr)


render("x:/project2016/woychik-ma-rnaseq/scripts/comparing_DEGs_v1.rmd", 
       output_file = "x:/project2016/woychik-ma-rnaseq/comparison/comparing_DEGs_v1.html")



render("x:/project2016/woychik-ma-rnaseq/scripts/comparing_DEGs_v2.rmd", 
       output_file = "x:/project2016/woychik-ma-rnaseq/comparison/comparing_DEGs_v1.html")

rmarkdown::render("x:/project2016/GoldenGateProject/RmdProj1/proj1_in_pdf.rmd", params = list (brfStatus = "breastfeeding",
                                                     sesStatus="social-status", 
                                                     genStatus="gender",
                                                     regStatus="region"), 
       output_file = "x:/project2016/GoldenGateProject/RmdProj1/proj1_in_pdf_02.pdf")

params$genStatus

setwd("X:/project2016/GSA_by_R")

render("GSA_4_Chen_v2.rmd", output_file = "GSA_4_Chen_v2_v1.html")


##  Let's see what is going on
render("GSA_4_Chen_data2_invest.rmd", output_file = "GSA_4_Chen_v2_original.html")

# This produce GSA analysis using the original 16 gene sets from David
render("GSA_4_Chen_data2_invest.rmd",  params = list(dbFile = "cleaningData/Knowledge_based_gene_sets_20140627.gmt"),
       output_file = "GSA_4_Chen_data1_orig.html")


# This produce GSA analysis using the original 16 gene sets from David but with slighly modification
# By adding "na" category
# Adding a new line at the end of gene set file
# Modifying the same Sox2xx to V2Sox2
render("GSA_4_Chen_data2_invest.rmd",  params = list(dbFile = "cleaningData/genset_1.txt"),
       output_file = "GSA_4_Chen_data1_mod.html")

# The results from two lists (dataset 1) seemed to be the same!


##  For second db sets by appending an extra set at the end


render("GSA_4_Chen_data2_invest.rmd",  params = list(dbFile = "cleaningData/Knowledge_based_gene_sets_20161005.gmt"),
       output_file = "GSA_4_Chen_data2_orig.html")
## The results are all messed up!!


## Manually appended the gene set on top of 20140627.gmt, Missing the last one
render("GSA_4_Chen_data2_invest.rmd",  params = list(dbFile = "cleaningData/temp.txt"),
       output_file = "GSA_4_Chen_data2_added_last_geneset.html")


## Manually modified temp.txt by creating a new line, NO problem
render("GSA_4_Chen_data2_invest.rmd", params = list(dbFile = "cleaningData/temp2.txt"),
       output_file = "GSA_4_Chen_v2_w_newline.html")



## Manually modified temp.txt by creating a new line, 
render("GSA_4_Chen_data2_invest.rmd", params = list(dbFile = "cleaningData/temp3.txt"),
       output_file = "GSA_4_Chen_v2_w_temp3.html")

render("GSA_4_Chen_data2_invest.rmd", output_file = "GSA_4_Chen_v2_cleaned_from_1.html")

render("GSA_4_Chen_data2_invest.rmd", output_file = "GSA_4_Chen_v2_cleaned_appended.html")
render("GSA_4_Chen_data2_invest.rmd", output_file = "GSA_4_Chen_v2_cleaned_appended2.html")

render("GSA_4_Chen_data2_invest.rmd",  params = list(dbFile = "cleaningData/temp4.txt") , 
       output_file = "GSA_4_Chen_v2_cleaned_appended3.html")


render("GSA_4_Chen_data2_invest.rmd", output_file = "GSA_4_Chen_v2_cleaned_appended4.html") 
## original dataset with "new line" added in the end

render("GSA_4_Chen_data2_invest.rmd",  params = list(dbFile = "cleaningData/temp4.txt") , 
       output_file = "GSA_4_Chen_v2_cleaned_appended4.html") 