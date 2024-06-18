library(tidyverse)
library(survival)
library(survminer)

## ---- Read in time to clad & censor data -----
KM = read.csv("../data/merged_KM_data.csv")
KM$microbiome_pid = as.character(KM$microbiome_pid)
KM = KM[,c("microbiome_pid", "clad_free_time", "status")]

meta = read.csv("../data/BAL_meta_filtered.csv", row.names = 1)

merged = inner_join(KM, meta, by = c("microbiome_pid" = "SubjectID")) %>% unique()
rownames(merged) = merged$SampleID

genus = read.csv("../data/BAL_genus_filtered.csv", row.names = 1)
genus = genus[ rownames(merged),]
percent = genus / rowSums(genus)

pseudomonas = percent[, "g__Pseudomonas"]

merged$pseudomonas = percent[, "g__Pseudomonas"] >= .1
merged$staph= percent[, "g__Staphylococcus"] >= .1

## ----- get on a per-patient basis ------
genera_by_patient = merged[, c("microbiome_pid", "clad_free_time", "status", "pseudomonas", "staph")] %>%
                    unique() %>%
                    group_by( microbiome_pid) %>% 
                    summarize( any_pseudomonas = any(pseudomonas), any_staph = any(staph)) %>%
                    ungroup() %>% as.data.frame()
KM_by_patient = left_join(genera_by_patient, merged[, c("microbiome_pid", "clad_free_time", "status")],
                            by = "microbiome_pid")

fit0 = survfit(Surv(clad_free_time, status)~1, data = KM_by_patient)
fit1 = survfit(Surv(clad_free_time, status)~any_pseudomonas, data = KM_by_patient)
fit2 = survfit(Surv(clad_free_time, status)~any_staph, data = KM_by_patient)

pdf("../figures/survival_psa_staph.pdf")
ggsurvplot(fit0, conf.int=TRUE, pval=TRUE)
ggsurvplot(fit1, conf.int=TRUE, pval=TRUE)
ggsurvplot(fit2, conf.int=TRUE, pval=TRUE)
dev.off()


coxph(Surv(clad_free_time, status) ~ any_pseudomonas, data = KM_by_patient)

## --------------- try all abundant taxa ---------------------
presence = genus > 5
presence = presence[ , colSums(presence) > ( ncol(presence) * .2 )]

cox_PH_list = lapply(colnames(presence), function(taxa){
    KM = merged
    KM$presence = presence[rownames(merged), taxa]
    by_patient = KM[, c("microbiome_pid", "clad_free_time", "status", "presence")] %>%
                    unique() %>%
                    group_by( microbiome_pid ) %>%
                    summarize(any_present = any(presence)) %>%
                    ungroup() %>% as.data.frame()
    KM_by_patient = left_join( by_patient, KM[, c("microbiome_pid", "clad_free_time", "status")],
                                by = "microbiome_pid") %>% unique()
    cox = summary(coxph(Surv(clad_free_time, status) ~ any_present, data = KM_by_patient))
    return(data.frame( taxa, coef = cox$coefficient[1], p = cox$coefficient[5]))
})
cox_PH = do.call("rbind", cox_PH_list)
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
    KM_by_patient = left_join( by_patient, KM[, c("microbiome_pid", "clad_free_time", "status")],
                                by = "microbiome_pid") %>% unique()
    cox = summary(coxph(Surv(clad_free_time, status) ~ any_present, data = KM_by_patient))
    return(data.frame( taxa, coef = cox$coefficient[1], p = cox$coefficient[5]))
})
cox_PH = do.call("rbind", cox_PH_list)
cox_PH$FDR = p.adjust(cox_PH$p, "BH")
arrange(cox_PH, p)



## ------------- plot ----------------
merged$hi_prevotella= percent[, "g__Prevotella"] >= .3
merged$hi_strep = percent[, "g__Streptococcus"] >= .3
genera_by_patient = merged[, c("microbiome_pid", "clad_free_time", "status",
                                                "hi_prevotella", "hi_strep")] %>%
                    unique() %>%
                    group_by( microbiome_pid) %>%
                    summarize(hi_prevotella = any(hi_prevotella), hi_strep = any(hi_strep)) %>%
                    ungroup() %>% as.data.frame()
KM_by_patient = left_join(genera_by_patient, merged[, c("microbiome_pid", "clad_free_time", "status")],
                            by = "microbiome_pid") %>% unique()
KM_by_patient$class = ifelse(KM_by_patient$hi_prevotella & !(KM_by_patient$hi_strep), "High_Prev.",
                       ifelse(KM_by_patient$hi_strep & !(KM_by_patient$hi_prevotella), "High_Strep.",
                        ifelse(KM_by_patient$hi_strep & KM_by_patient$hi_prevotella, "Both", "Neither")))
fit1 = survfit(Surv(clad_free_time, status)~class, data = KM_by_patient)
fit2 = survfit(Surv(clad_free_time, status)~hi_strep, data = KM_by_patient)
fit3 = survfit(Surv(clad_free_time, status)~hi_prevotella, data = KM_by_patient)

pdf("../figures/high_strep_prevotella.pdf")
ggsurvplot(fit1, conf.int=TRUE, pval=TRUE)
ggsurvplot(fit2, conf.int=TRUE, pval=TRUE)
ggsurvplot(fit3, conf.int=TRUE, pval=TRUE)
dev.off()
summary(coxph( Surv(clad_free_time, status) ~ class, data = KM_by_patient))

cox_PH$FDR = p.adjust(cox_PH$p, "BH")
arrange(cox_PH, p)

## --------------- highly abundant taxa ---------------------
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
    KM_by_patient = left_join( by_patient, KM[, c("microbiome_pid", "clad_free_time", "status")],
                                by = "microbiome_pid") %>% unique()
    cox = summary(coxph(Surv(clad_free_time, status) ~ any_present, data = KM_by_patient))
    return(data.frame( taxa, coef = cox$coefficient[1], p = cox$coefficient[5]))
})
cox_PH = do.call("rbind", cox_PH_list)
cox_PH$FDR = p.adjust(cox_PH$p, "BH")
arrange(cox_PH, p)


##  plot prevotella
merged$hi_prevotella= percent[, "g__Prevotella"] >= .3
merged$hi_strep = percent[, "g__Streptococcus"] >= .3
genera_by_patient = merged[, c("microbiome_pid", "clad_free_time", "status", 
                                                "hi_prevotella", "hi_strep")] %>%
                    unique() %>%
                    group_by( microbiome_pid) %>% 
                    summarize(hi_prevotella = any(hi_prevotella), hi_strep = any(hi_strep)) %>%
                    ungroup() %>% as.data.frame()
KM_by_patient = left_join(genera_by_patient, merged[, c("microbiome_pid", "clad_free_time", "status")],
                            by = "microbiome_pid") %>% unique()
KM_by_patient$class = ifelse(KM_by_patient$hi_prevotella & !(KM_by_patient$hi_strep), "High_Prev.",
                       ifelse(KM_by_patient$hi_strep & !(KM_by_patient$hi_prevotella), "High_Strep.",
                        ifelse(KM_by_patient$hi_strep & KM_by_patient$hi_prevotella, "Both", "Neither")))
fit1 = survfit(Surv(clad_free_time, status)~class, data = KM_by_patient)
fit2 = survfit(Surv(clad_free_time, status)~hi_strep, data = KM_by_patient)
fit3 = survfit(Surv(clad_free_time, status)~hi_prevotella, data = KM_by_patient)

pdf("../figures/high_strep.pdf")
ggsurvplot(fit2, conf.int=TRUE, pval=TRUE, break.x.by = 365,
            palette = c("#0a6102", "#ffa217"), xlim = c(0, 2922))
dev.off()
pdf("../figures/high_prevotella.pdf")
ggsurvplot(fit3, conf.int=TRUE, pval=TRUE, break.x.by = 365,
            palette = c("#0a6102", "#ffa217"), xlim = c(0, 2922))
dev.off()

summary(coxph( Surv(clad_free_time, status) ~ class, data = KM_by_patient))

###---------------------------------------------------------------
merged$prev_strep = (.001 + percent[,"g__Prevotella"]) / (.001 + percent[,"g__Streptococcus"] )

prev_strep_ph = lapply(unique(merged$timepoint), function(tp){

    filt = merged[ merged$timepoint == tp,]
    cox = summary(coxph(Surv(clad_free_time, status) ~ prev_strep, data = filt))
    return(data.frame( tp, coef = cox$coefficient[1], p = cox$coefficient[5]))
})
prev_strep_ph = do.call('rbind', prev_strep_ph)
prev_strep_ph

M6 = merged[ merged$timepoint == "6M",]
M6$class = cut(M6$prev_strep, quantile(M6$prev_strep))
summary(coxph(Surv(clad_free_time, status) ~ class, data = M6))

sub = sample_n(M6, 100)
summary(coxph(Surv(clad_free_time, status) ~ prev_strep, data = sub))
summary(coxph(Surv(clad_free_time, status) ~ prev_strep, data = filter(M6, prev_strep < 100)))


library(ROCR)
pred = prediction(M6$prev_strep, M6$status)
perf = performance(pred, "tpr", "fpr")
pdf("../figures/tmp_roc.pdf")
plot(perf)
dev.off()

