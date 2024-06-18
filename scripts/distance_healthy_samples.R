library(tidyverse)
library(vegan)
library(survival)

## get the UniFrac distance between each timepont and the last sample from the same subject
## Then get total distance / total time

meta = read.csv("../data/BAL_meta_filtered.csv", row.names = 1)
meta$timepoint = factor(meta$timepoint, levels = c("H", "D", "R", "DC", "6W", "3M", "6M", "12M"))
weighted = readRDS("../data/weighed_unifrac_no_contam.rds")

## only get the overlaps - the unifract distance will drop some more samples because
## we removed the contaminants and now more are <1000 reads
meta = meta[ meta$SampleID %in% rownames(weighted),]


## ----- Distance to healthy control BAL samples  ---------------------------
healthy = meta[ meta$timepoint == "H" & !(is.na(meta$timepoint)), "SampleID"]
transplant = meta[ meta$timepoint != "H" & !(is.na(meta$timepoint)), "SampleID"]

distance_to_healthy = rowMeans(weighted[ transplant, healthy ])

transplant_meta = meta[transplant,]
transplant_meta$distance_to_healthy = distance_to_healthy[ rownames(transplant_meta)]
wilcox.test(distance_to_healthy ~ timepoint, data = filter(transplant_meta, timepoint %in% c("D", "R")))

pdf("../figures/distance_to_healthy.pdf", width = 12, height =8)
ggplot(transplant_meta, aes(x = timepoint, y = distance_to_healthy) ) +
    theme_classic() +
    geom_point() +
    geom_line(aes(group = SubjectID), linewidth = .5) +
    geom_boxplot(outlier.alpha = 0, alpha = .5)
dev.off()

## ----- Does distance to healthy have any effect on outcome? ----------------
KM_data = read.csv("../data/merged_KM_data.csv")

dist = data.frame(sample = names(distance_to_healthy), distance_to_healthy,
                    quantile_dist = cut(distance_to_healthy, breaks = 4,
                                        labels = c("low", "mid_low", "mid_hi", "hi")))
dist = merge( dist, meta, by.x = "sample", by.y = "SampleID")
KM_dist = merge(dist, KM_data, by.x = "SubjectID", by.y = "microbiome_pid")

pvals = sapply(unique(dist$timepoint), function(tp) {
        filt = KM_dist %>% filter(timepoint == tp)
        print(tp)
        cox = coxph(Surv(clad_free_time, status) ~ distance_to_healthy + age, data = filt) 
        print(cox)
        return(summary(cox)$coefficients[1, 5])})







