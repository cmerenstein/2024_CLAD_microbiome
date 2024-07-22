library(tidyverse)
library(survival)
library(survminer)

## ---- Read in time to clad & censor data -----
KM = read.csv("../data/merged_KM_data.csv")
KM = KM[,c("microbiome_pid", "age", "status", "clad_free_time")]
KM$microbiome_pid = as.character(KM$microbiome_pid)

meta = read.csv("../data/BAL_meta_filtered.csv", row.names = 1)

merged = inner_join(KM, meta, by = c("microbiome_pid" = "SubjectID")) %>% unique()
rownames(merged) = merged$SampleID

genus = read.csv("../data/BAL_genus_filtered.csv", row.names = 1)
genus = genus[ rownames(merged),]
percent = genus / rowSums(genus)

## ------ specifically check pseudomonas ------------
merged$pseudomonas = percent[ , "g__Pseudomonas"]
genera_by_patient = merged[, c("microbiome_pid", "clad_free_time", "status", "pseudomonas")] %>%
                    unique() %>%
                    group_by( microbiome_pid) %>% 
                    summarize(hi_pseudomonas = any(pseudomonas > .3)) %>%
                    ungroup() %>% as.data.frame()

KM_by_patient = left_join(genera_by_patient,
                             merged[,c("microbiome_pid", "age", "status", "clad_free_time")], 
                            by = "microbiome_pid") %>% unique()

fit2 = survfit(Surv(clad_free_time, status)~hi_pseudomonas, data = KM_by_patient)

pdf("../figures/high_pseudomonas.pdf")
ggsurvplot(fit2, conf.int=TRUE, pval=TRUE, break.x.by = 365,
            palette = c("#0a6102", "#ffa217"), xlim = c(0, 2922))
dev.off()

summary(coxph( Surv(clad_free_time, status) ~ hi_pseudomonas+ age, data = KM_by_patient))

###---------------------------------------------------------------
pseudomonas_ph = lapply(unique(merged$timepoint), function(tp){

    filt = merged[ merged$timepoint == tp,]
    cox = summary(coxph(Surv(clad_free_time, status) ~ pseudomonas + age, data = filt))
    print(cox)
    return(data.frame( tp, HR = cox$coefficient[1, 2], p = cox$coefficient[1, 5]))
})
pseudomonas_ph = do.call('rbind', pseudomonas_ph)
pseudomonas_ph

lapply(unique(merged$timepoint), function(tp){

    filt = merged[ merged$timepoint == tp,]
    wx = wilcox.test( filt$pseudomonas ~ filt$status)
    print(wx)
})





