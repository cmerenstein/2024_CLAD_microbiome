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

## ------------- plot ----------------
merged$hi_prevotella= percent[, "g__Prevotella"] >= .3
merged$hi_strep = percent[, "g__Streptococcus"] >= .3
genera_by_patient = merged[, c("microbiome_pid", "clad_free_time", "status", 
                                                "hi_prevotella", "hi_strep")] %>%
                    unique() %>%
                    group_by( microbiome_pid) %>% 
                    summarize(hi_prevotella = any(hi_prevotella), hi_strep = any(hi_strep)) %>%
                    ungroup() %>% as.data.frame()

KM_by_patient = left_join(genera_by_patient,
                             merged[,c("microbiome_pid", "age", "status", "clad_free_time")], 
                            by = "microbiome_pid") %>% unique()

KM_by_patient$class = ifelse(KM_by_patient$hi_prevotella & !(KM_by_patient$hi_strep), "High_Prev.",
                       ifelse(KM_by_patient$hi_strep & !(KM_by_patient$hi_prevotella), "High_Strep.",
                        ifelse(KM_by_patient$hi_strep & KM_by_patient$hi_prevotella, "Both", "Neither")))

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
summary(coxph( Surv(clad_free_time, status) ~ hi_strep, data = KM_by_patient))
summary(coxph( Surv(clad_free_time, status) ~ hi_prevotella, data = KM_by_patient))
summary(coxph( Surv(clad_free_time, status) ~ hi_strep + age, data = KM_by_patient))
summary(coxph( Surv(clad_free_time, status) ~ hi_prevotella + age, data = KM_by_patient))

###---------------------------------------------------------------
merged$prev_strep = log10( 1 + genus[,"g__Streptococcus"]) - 
                    log10(1 + genus[,"g__Prevotella"])

prev_strep_ph = lapply(unique(merged$timepoint), function(tp){

    filt = merged[ merged$timepoint == tp,]
    cox = summary(coxph(Surv(clad_free_time, status) ~ prev_strep, data = filt))
    return(data.frame( tp, HR = cox$coefficient[2], p = cox$coefficient[5]))
})
prev_strep_ph = do.call('rbind', prev_strep_ph)
prev_strep_ph

W6 = merged[ merged$timepoint == "6W",]
W6$class = cut(W6$prev_strep, quantile(W6$prev_strep, c(0,.33, .67,1)))
summary(coxph(Surv(clad_free_time, status) ~ class, data = W6))
summary(coxph(Surv(clad_free_time, status) ~ prev_strep + age, data = W6))
wilcox.test(prev_strep ~ status, data = W6)

pdf("../figures/week_6_prev_strep.pdf")
ggsurvplot( survfit(Surv(clad_free_time, status) ~ class, data = W6), 
            xlim = c(0, 2922), conf.int=TRUE, pval=TRUE, break.x.by = 365)
dev.off()

## make sure finding is robust to outlyers
shuffle = sapply(seq(1:100), function(i){
    sub = sample_n(W6, 100)
    cox = summary(coxph(Surv(clad_free_time, status) ~ prev_strep, data = sub))
    return(cox$coefficient[5])
    })
sum(shuffle < .1)

sapply(seq(1:132), function(i){
    sub = W6[-i,]
    cox = summary(coxph(Surv(clad_free_time, status) ~ prev_strep, data = sub))
    return(cox$coefficient[5])
    })


###-----------------------------------------------------------
merged$time_from_transplant = as.numeric( as.Date(merged$DateCollected, "%m/%d/%Y") -
                                as.Date(merged$transplant_date, "%m/%d/%Y"))
merged$CLAD_status = c("non-CLAD", "CLAD")[1 + merged$status]

pdf("../figures/prev_strep_ratio_time_from_transplant.pdf")
ggplot(filter(merged, time_from_transplant > -2), 
        aes(x = time_from_transplant, y = prev_strep, color = CLAD_status)) +
    theme_classic() +
    geom_point() +
#    geom_line(aes(group = microbiome_pid), alpha = 0.2, linewidth = .5) +
    scale_color_manual(values = c("#f75965", "#a2d6fa")) +
    stat_smooth() +
    xlab("Time from transplant") +
    ylab("Streptococcus / Prevotella log ratio") +
    guides(color = guide_legend(title = "CLAD status")) +
    theme(text = element_text(size = 16)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    xlim(c(0, 400)) + 
    theme(legend.position = c(0.8, 0.2))
dev.off()







