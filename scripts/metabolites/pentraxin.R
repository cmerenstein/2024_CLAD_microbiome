library(tidyverse)
library(survival)
library(survminer)


KM = read.csv("../../data/merged_KM_data.csv")
meta = read.csv("../../data/BAL_meta_unfiltered.csv", row.names = 1)
meta$alt_id = sapply(meta$SampleID, function(d) {strsplit(d, "\\.")[[1]][1]})

merged = merge(KM, meta, by.x = "microbiome_pid", by.y = "SubjectID")

pentraxin = read.csv("../../data/metabolites/pentraxin_clean.csv")[1:4]
pentraxin$Analyte = "Pentraxin"
colnames(pentraxin)[2] = "Concentration"
pentraxin = pentraxin[ !(is.na(pentraxin$Concentration)),]

pentraxin = merge(pentraxin, merged, by.x = "sample", by.y = "alt_id")
summary(coxph(Surv(clad_free_time, status)~Concentration + age, data = pentraxin))

