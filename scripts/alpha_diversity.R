library(vegan)
library(ape)
library(tidyverse)
library(survival)
library(survminer)


## ---------- Read in all metadata 
BAL = read.csv("../data/BAL_meta_filtered.csv", row.names = 1)
genus = read.csv("../data/genus_no_contam.csv", row.names = 1)
genus = genus[rowSums(genus) >= 1000,]

percent = genus / rowSums(genus)

simpson = diversity(percent, "simpson")
shannon = diversity(percent, "shannon")
chao = estimateR(genus)["S.chao1",]

diversity = rbind( data.frame(sample = names(simpson), value = simpson, metric = "simpson"),
                    data.frame(sample = names(shannon), value = shannon, metric = "shannon"),
                    data.frame(sample = names(chao), value = chao, metric = "chao"))

#write.csv(diversity, "../data/alpha_diversity.csv")

## ------ QC by batch and depth -----------------------
diversity_meta = merge(diversity, BAL, by.x = "sample", by.y = "SampleID")

pdf("../figures/diversity_QC.pdf")

ggplot(diversity_meta, aes(x = as.factor(run), y = value) ) + 
    facet_wrap(~metric, scale = "free") + 
    theme_classic() + 
    geom_violin()

ggplot(diversity_meta, aes(x = depth, y = value)) + 
    facet_wrap(~metric, scale = "free") + 
    theme_classic() + 
    geom_point() +
    scale_x_log10()
dev.off()


## ------ Diversity by timepoint -----
KM = read.csv("../data/merged_KM_data.csv")
KM$microbiome_pid = as.character(KM$microbiome_pid)

transplant = diversity_meta[!(is.na(diversity_meta$timepoint)) & diversity_meta$timepoint != "H",]
transplant$timepoint = factor(transplant$timepoint,
                                levels = c("H", "D", "R", "DC", "6W", "3M", "6M", "12M"))
transplant$run = as.factor(transplant$run)
merged = merge(KM, transplant, by.x = "microbiome_pid", by.y = "SubjectID")

cox_PH_timepoint = lapply(unique(transplant$timepoint), function(tp){
                        print(tp)
                        filt = filter(merged, timepoint == tp & metric == "simpson")
                        cox = coxph(Surv(clad_free_time, status) ~ value + age, data = filt)
                        return(data.frame(tp, p = summary(cox)$coefficients[1,5], 
                                          OR = summary(cox)$coefficients[1,2])) })
do.call("rbind", cox_PH_timepoint)


pdf("../figures/diversity_CLAD.pdf")

for( tp in c("D", "R", "DC", "6W", "3M", "6M", "12M")){
    filtered = filter(merged, timepoint == tp & metric == "simpson")
    quantiles = c(-Inf, quantile(filtered$value, c(.25, .75)), Inf)
    filtered$diversity_bins = cut(filtered$value, quantiles, labels = c("low", "mid", "high"))
    fit = survfit(Surv(clad_free_time, status) ~ diversity_bins, data = filtered)
    print(ggsurvplot(fit, conf.int = T, pval = T, 
            break.x.by = 365, xlim = c(0, 2922), title = tp))
}
dev.off()


















