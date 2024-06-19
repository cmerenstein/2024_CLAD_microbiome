library(tidyverse)
library(survival)
library(survminer)

## ---- Read in time to clad & censor data -----
KM = read.csv("../data/merged_KM_data.csv")
KM$microbiome_pid = as.character(KM$microbiome_pid)

meta = read.csv("../data/BAL_meta_filtered.csv", row.names = 1)

merged = inner_join(KM, meta, by = c("microbiome_pid" = "SubjectID")) %>% unique()
rownames(merged) = merged$SampleID

genus = read.csv("../data/BAL_genus_filtered.csv", row.names = 1)
genus = genus[ rownames(merged),]
percent = genus / rowSums(genus)

## --------------- try all abundant taxa ---------------------
presence = genus > 5 ## threshold of 5 reads
presence = presence[ , colSums(presence) > ( ncol(presence) * .2 )] ## any taxa in at least 20% of samples

cox_PH_list = lapply(colnames(presence), function(taxa){
    KM = merged
    KM$presence = presence[rownames(merged), taxa]
    by_patient = KM[, c("microbiome_pid", "clad_free_time", "status", "presence")] %>%
                    unique() %>%
                    group_by( microbiome_pid ) %>%
                    summarize(any_present = any(presence)) %>%
                    ungroup() %>% as.data.frame()
    KM_by_patient = left_join( by_patient, KM[, c("microbiome_pid", "clad_free_time","age", "status")],
                                by = "microbiome_pid") %>% unique()
    cox = summary(coxph(Surv(clad_free_time, status) ~ any_present + age, data = KM_by_patient))
    return(data.frame( taxa, coef = cox$coefficient[1], p = cox$coefficient[5]))
})
cox_PH = do.call("rbind", cox_PH_list)
cox_PH$fdr = p.adjust(cox_PH$p, "BH")
arrange(cox_PH, p)

## ------------------- Just highly abundant -----------------------
abundant = percent > .3
abundant = abundant[ , colSums(abundant) > ( ncol(abundant) * .1 )]

cox_PH_list = lapply(colnames(abundant), function(taxa){
    KM = merged
    KM$abundant = abundant[rownames(merged), taxa]
    by_patient = KM[, c("microbiome_pid", "clad_free_time", "status", "abundant")] %>%
                    unique() %>%
                    group_by( microbiome_pid ) %>%
                    summarize(any_present = any(abundant)) %>%
                    ungroup() %>% as.data.frame()
    KM_by_patient = left_join( by_patient, KM[, c("microbiome_pid", "clad_free_time","age", "status")],
                                by = "microbiome_pid") %>% unique()
    cox = summary(coxph(Surv(clad_free_time, status) ~ any_present + age, data = KM_by_patient))
    return(data.frame( taxa, coef = cox$coefficient[1], p = cox$coefficient[5]))
})
cox_PH = do.call("rbind", cox_PH_list)
cox_PH$FDR = p.adjust(cox_PH$p, "BH")
arrange(cox_PH, p)

