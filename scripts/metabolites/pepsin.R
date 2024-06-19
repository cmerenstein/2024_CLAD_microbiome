library(tidyverse)
library(survival)
library(survminer)


KM = read.csv("../../data/merged_KM_data.csv")
meta = read.csv("../../data/BAL_meta_unfiltered.csv", row.names = 1)
meta$alt_id = sapply(meta$SampleID, function(d) {strsplit(d, "\\.")[[1]][1]})

merged = merge(KM, meta, by.x = "microbiome_pid", by.y = "SubjectID")

pepsin = read.csv("../../data/metabolites/pepsin_clean.csv")[1:4]
pepsin$Analyte = "Pepsin"
colnames(pepsin)[2] = "Concentration"
colnames(pepsin)[4] = "Std.Dev"

pepsin = merge(pepsin, merged, by.x = "sample", by.y = "alt_id")
summary(coxph(Surv(clad_free_time, status)~Concentration + age, data = pepsin))

pepsin$quantile = c("low", "mid", "high")[as.numeric(cut(pepsin$Concentration, 3))]
pepsin$presence = pepsin$Concentration > 0
fit1 = survfit(Surv(clad_free_time, status)~presence, data = pepsin)

pdf("../../figures/luminex/pepsin_survival.pdf")
ggsurvplot(fit1, conf.int=TRUE, pval=TRUE)
dev.off()


