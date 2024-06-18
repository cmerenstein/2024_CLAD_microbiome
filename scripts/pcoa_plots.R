library(tidyverse)
library(vegan)
library(ape)
library(survival)
library(survminer)

## get metadata and unifrac
meta = read.csv("../data/BAL_meta_filtered.csv", row.names = 1)
meta = meta[ !(is.na(meta$timepoint)),]
meta$Time = as.Date(meta$DateCollected, format = "%m/%d/%Y")

## merge R and RT
meta$timepoint[meta$timepoint == "RT"] <- "R"
meta$timepoint[meta$timepoint == "H"] <- "Healthy"

## read unifrac distances
rare_weighted = readRDS("../data/weighed_unifrac_no_contam.rds")

## only get the overlaps - the unifract distance will drop some more samples because
## we removed the contaminants and now more are <1000 reads
meta = meta[ meta$SampleID %in% rownames(rare_weighted),]
rare_weighted = rare_weighted[meta$SampleID, meta$SampleID]

## Clad onset data
KM = read.csv("../data/merged_KM_data.csv", row.names = 1)
CLAD_status = KM[,c("microbiome_pid", "status", "clad_free_time")]

meta_CLAD = merge(meta, CLAD_status, by.x = "SubjectID", by.y = "microbiome_pid", all.x = T)
rownames(meta_CLAD) = meta_CLAD$SampleID

## ------- Generate PCoA ordinates ---------------------
set.seed(19002)
ordinates = pcoa(rare_weighted)

str(ordinates) ## look at relative correlation for variance explained
vectors = ordinates$vectors

pcoa_df = data.frame(sample = rownames(vectors), Axis.1 = vectors[,1],
                    Axis.2 = vectors[,2], Axis.3 = vectors[,3], Axis.4 = vectors[,4],
                    run = as.factor(meta_CLAD[rownames(vectors), "run"]),
                    CLAD = as.factor(meta_CLAD[rownames(vectors), "status"]),
                    timepoint = factor(meta_CLAD[rownames(vectors), "timepoint"],
                        levels = c("Healthy", "D", "R", "DC", "6W", "3M", "6M", "12M")))
pcoa_df_filter = pcoa_df[ !(is.na(pcoa_df$CLAD)) | pcoa_df$timepoint == "Healthy",]

pdf("../figures/figure_1_weighted_unifrac_pcoa.pdf")
ggplot(pcoa_df_filter, aes(x = Axis.1, y = Axis.2, color = CLAD)) + 
    geom_point() + 
    stat_ellipse() +
    theme_classic() +
    theme(text = element_text(size = 14)) + 
    facet_wrap(~timepoint, ncol = 3) +
    scale_color_manual(values = c("#a2d6fa", "#f75965"))
dev.off() 

## --------- permanova -----------------------------

transplant = meta_CLAD %>% filter(timepoint != "Healthy")
null = lapply(unique(transplant$timepoint), function(tp){

    filt = transplant %>% filter(timepoint == tp) %>%
                            filter(!(is.na(status)))
    filt_dist = as.matrix(rare_weighted)[rownames(filt), rownames(filt)]
    print(tp)
    print(adonis(filt_dist~ status, data = filt)$aov.tab)
})

## --------- Cox PH ----------------------------
cox_list = lapply(unique(transplant$timepoint), function(tp){
    filt = transplant %>% filter(timepoint == tp) %>%
                            filter(!(is.na(status)))
    filt$PC1 = vectors[rownames(filt), "Axis.1"] 
    filt$PC2 = vectors[rownames(filt), "Axis.2"] 
    filt$PC3 = vectors[rownames(filt), "Axis.3"] 
    filt$PC4 = vectors[rownames(filt), "Axis.4"] 
    cox = summary(coxph(Surv(clad_free_time, status)~ PC1 + PC2 + PC3 + PC4, data = filt))    
    df = rbind( data.frame( tp, axis = "PC1", p = cox$coefficients[1,5], coef = cox$coefficients[1,2]),
           data.frame( tp, axis = "PC2", p = cox$coefficients[2,5], coef = cox$coefficients[2,2]),
           data.frame( tp, axis = "PC3", p = cox$coefficients[3,5], coef = cox$coefficients[3,2]),
           data.frame( tp, axis = "PC4", p = cox$coefficients[4,5], coef = cox$coefficients[4,2]))
    return(df)
})
cox_df = do.call("rbind", cox_list)
cox_df$fdr = p.adjust(cox_df$p, "BH")
arrange(cox_df, p)

## ------------- See what's correlated to what --------------------
genus = read.csv("../data/BAL_genus_filtered.csv", row.names = 1)
genus = genus[ rownames(vectors),]
percent = genus / rowSums(genus)
common = percent[, colMeans(percent) > .05]

fit = envfit(vectors, common)
genus_fit = fit$vectors$arrows %>% as.data.frame
genus_fit$genus = rownames(genus_fit)

pdf("../figures/pcoa_loadings.pdf", width = 4, height = 4)
ggplot(pcoa_df_filter, aes(x = Axis.1, y = Axis.2)) +
    theme_classic() + 
    geom_point(alpha = 0) +
    geom_segment(data = genus_fit, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
                    arrow = arrow(length = unit(0.2, "cm")) ) +
    geom_text(data = genus_fit, aes(x = Axis.1 / 2, y = Axis.2 / 2, label = genus))
dev.off()




