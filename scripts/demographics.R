library(tidyverse)
library(survival)
library(survminer)

## ---- Read in time to clad & censor data -----
KM = read.csv("../data/merged_KM_data.csv")
KM$microbiome_pid = as.character(KM$microbiome_pid)
KM = KM[,c("microbiome_pid", "clad_free_time", "status", "days_from_transplant")]

## ---- Read in metadata and merge ----------------
meta = read.csv("../data/BAL_meta_filtered.csv", row.names = 1)
merged = inner_join(KM, meta, by = c("microbiome_pid" = "SubjectID")) %>% unique()
rownames(merged) = merged$SampleID

## --- Demographics -----------------------------
demo = read.csv("../data/2024-05-03additionalmicrodata.csv")
demo_CLAD = merge(KM, demo, by = "microbiome_pid")
demo_CLAD$age_bin = cut(demo_CLAD$age, 4)
demo_CLAD$BMI_bin = cut(demo_CLAD$BMI, 4)

fit0 = survfit(Surv(clad_free_time, status)~1, data = demo_CLAD)
fit1 = survfit(Surv(clad_free_time, status)~gender, data = demo_CLAD)
fit2 = survfit(Surv(clad_free_time, status)~age_bin, data = demo_CLAD)
fit3 = survfit(Surv(clad_free_time, status)~BMI_bin, data = demo_CLAD)

pdf("../figures/demographics_survival.pdf")
ggsurvplot(fit1, conf.int=TRUE, pval=TRUE, break.x.by = 365,
            palette = c("#0a6102", "#ffa217"), xlim = c(0, 2922))
ggsurvplot(fit2, conf.int=TRUE, pval=TRUE, break.x.by = 365,
             xlim = c(0, 2922))
ggsurvplot(fit3, conf.int=TRUE, pval=TRUE, break.x.by = 365,
             xlim = c(0, 2922))
dev.off()

coxph( Surv(clad_free_time, status) ~ gender, data = demo_CLAD)
coxph( Surv(clad_free_time, status) ~ age, data = demo_CLAD)
coxph( Surv(clad_free_time, status) ~ BMI, data = demo_CLAD)

## ---- means for table --------------------
group_by(demo_CLAD, status) %>% summarize(age = mean(age), BMI = mean(BMI), 
                followup = mean(days_from_transplant))
sd(demo_CLAD[ demo_CLAD$status == 1, "age"])
sd(demo_CLAD[ demo_CLAD$status == 0, "age"])
sd(demo_CLAD[ demo_CLAD$status == 1, "BMI"])
sd(demo_CLAD[ demo_CLAD$status == 0, "BMI"])
sd(demo_CLAD[ demo_CLAD$status == 1, "days_from_transplant"])
sd(demo_CLAD[ demo_CLAD$status == 0, "days_from_transplant"])

disease = unique(merged[,c("microbiome_pid", "status", "diagnosis", "transplant")])
table(disease[,c("diagnosis", "status")])




