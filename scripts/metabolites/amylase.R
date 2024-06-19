library(tidyverse)
library(survival)
library(survminer)


KM = read.csv("../../data/merged_KM_data.csv")
meta = read.csv("../../data/BAL_meta_unfiltered.csv", row.names = 1)
meta$alt_id = sapply(meta$SampleID, function(d) {strsplit(d, "\\.")[[1]][1]})

merged = merge(KM, meta, by.x = "microbiome_pid", by.y = "SubjectID")

amylase = read.csv("../../data/metabolites/amylase_clean.csv")
amylase$Analyte = "Amylase"
amylase = amylase[,c("sample", "Amylase.corrected", "Analyte")]

amylase = merge(amylase, merged, by.x = "sample", by.y = "alt_id")
summary(coxph(Surv(clad_free_time, status)~Amylase.corrected + age, data = amylase))

