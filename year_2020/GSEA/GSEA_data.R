overall.gender.diff.dt <- zoom.in.dt

colnames(overall.gender.diff.dt) <- sub (".mm9.genes.fpkm_tracking.fpkm.txt", "", colnames(overall.gender.diff.dt))

rownames(overall.gender.diff.dt) <-  getGeneList(rownames(overall.gender.diff.dt))
str(overall.gender.diff.dt)
file.2.save <- paste ("x:/project2013/", "DEG_overall_gender_diff.txt", sep = "")
write.table (overall.gender.diff.dt, file.2.save, sep = "\t",  row.names= TRUE, col.names=NA)

