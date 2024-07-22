library(tidyverse)
library(survival)
library(survminer)

## ---- Read in time to clad & censor data -----
KM = read.csv("../data/merged_KM_data.csv")
KM$microbiome_pid = as.character(KM$microbiome_pid)

meta = read.csv("../data/BAL_meta_filtered.csv", row.names = 1)

KM_meta = inner_join(KM, meta, by = c("microbiome_pid" = "SubjectID")) %>% unique()


## ----- read in PGD data -------------------
PGD = read.csv("../data/PGD_grades.csv")
PGD$max_PGD = apply(PGD[, c("pgd_day_0", "pgd_day_1", "pgd_day_2", "pgd_day_3")], 1, max, na.rm = T)

PGD$PGD_score = apply(PGD[, c("pgd_day_0", "pgd_day_1", "pgd_day_2", "pgd_day_3")], 1, function(scores){
                scores = scores[!(is.na(scores)) & !( scores == -99)]
                if( all(scores == 3)) { return( "PGD-3")}
                if( all(scores == 0)) { return( "PGD-0")}
                else { return("intermediate")}
            })

PGD_filt = PGD[, c("microbiome_pid", "max_PGD", "PGD_score")]
PGD_filt$PGD_score = factor(PGD_filt$PGD_score, levels = c("PGD-0", "intermediate", "PGD-3"))
write.csv(PGD_filt, "../data/composite_PGD_score.csv", row.names = F)

merged = merge(KM_meta, PGD_filt, by = "microbiome_pid")
rownames(merged) = merged$SampleID

genus = read.csv("../data/BAL_genus_filtered.csv", row.names = 1)
genus = genus[ rownames(merged),]
percent = genus / rowSums(genus)

merged$prev_strep = log10( 1 + genus[,"g__Prevotella"]) -
                    log10(1 + genus[,"g__Streptococcus"])

## ------- association and mediation ----------------------------
W6 = merged[ merged$timepoint == "6W",]
summary(coxph(Surv(clad_free_time, status) ~ prev_strep + age, data = W6))
summary(coxph(Surv(clad_free_time, status) ~ PGD_score + age, data = W6))
summary(coxph(Surv(clad_free_time, status) ~ prev_strep + PGD_score + age, data = W6))

clinfun::jonckheere.test(W6$prev_strep, as.numeric(W6$PGD_score), "increasing")
clinfun::jonckheere.test(merged[ merged$timepoint == "R", "prev_strep"], 
                         as.numeric(merged[ merged$timepoint == "R", "PGD_score"]), "increasing")

pdf("../figures/PGD_prev_strep_recipient.pdf")
ggplot( filter(merged, timepoint == "R"), aes(x = PGD_score, y = prev_strep)) + 
    theme_bw() +
    geom_boxplot() + 
    geom_jitter(height = 0, width = 0.33)
dev.off()



