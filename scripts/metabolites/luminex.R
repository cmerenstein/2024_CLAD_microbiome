library(tidyverse)
library(survival)
library(survminer)
library(vegan)

meta = read.csv("../../data/BAL_meta_unfiltered.csv", row.names = 1)
meta$alt_id = sapply(meta$SampleID, function(d) {strsplit(d, "\\.")[[1]][1]})

## ---------------
KM = read.csv("../../data/merged_KM_data.csv")

merged = merge(KM, meta, by.x = "microbiome_pid", by.y = "SubjectID") %>%
            select("microbiome_pid", "alt_id", "status", "clad_free_time", "age", "SampleID")

luminex = read.csv("../../data/metabolites/luminex_clean.csv")[, 1:5]
luminex_KM = left_join(luminex, merged, by = c("sample" = "alt_id"))

## -------- luminex KM list ------------
met_KM_list = lapply(unique(luminex_KM$Analyte), function(a){
            filt = filter(luminex_KM, Analyte == a)
            cox = summary(coxph(Surv(clad_free_time, status) ~ Concentration + age, data = filt))
            return(data.frame( metabolte = a, HR= cox$coefficients[1,2], p = cox$coefficients[1,5]))})
met_KM = do.call("rbind", met_KM_list)
met_KM$FDR = p.adjust(met_KM$p, method = "BH")

write.csv(met_KM, "coxPH_luminex.csv", row.names = F)
arrange(met_KM, p)

pdf("../../figures/IP10.pdf")
IP10 = luminex_KM[luminex_KM$Analyte == "IP-10" & !(is.na(luminex_KM$SampleID)),]
IP10$class = ifelse(IP10$Concentration > 710, "high", "normal") %>% factor(levels = c("normal", "high"))
ggsurvplot( survfit(Surv(clad_free_time, status) ~ class, data = IP10),
            conf.int = T, pval = T, break.x.by = 365,
            palette = c( "#ffa217", "#0a6102"), xlim = c(0, 2922))
dev.off()

summary(coxph( Surv(clad_free_time, status) ~ class + age, data = IP10))

## make sure it's not being driven by just a few samples
shuffle = sapply(seq(1:1000), function(i){
    sub = sample_n(IP10, 80)
    cox = summary(coxph(Surv(clad_free_time, status) ~ Concentration + age, data = sub))
    p = cox$coefficient[5]
    if( p > .1){ print( IP10[ !(IP10$sample %in% sub$sample),
                       c("sample", "Concentration", "clad_free_time")])}
    return(p)
    })
sum(shuffle < .1)



pdf("../../figures/null_survival.pdf")
ggsurvplot( survfit( Surv(clad_free_time, status) ~ 1, data = KM),
            palette = c("black"), break.x.by = (365 * 2), 
            xlim = c(0, 3650), surv.median.line = "hv")
dev.off()


## --------------------------------------------------------------
luminex_long = spread(luminex[, c("sample", "Analyte", "Concentration")], 
                            "Analyte", "Concentration")
luminex_mat = as.matrix(luminex_long[, 2:ncol(luminex_long)])
rownames(luminex_mat) = luminex_long$sample
luminex_pca = prcomp(luminex_mat)

CLAD_status = KM[,c("microbiome_pid", "status")]
plot_pca = data.frame(sample = rownames(luminex_pca$x), luminex_pca$x[,1:5]) %>%
                left_join(meta, by = c("sample" = "alt_id"))
plot_pca = merge(plot_pca, CLAD_status, by.x = "SubjectID", by.y = "microbiome_pid",
                 all.x = T, all.y = F) %>%
                 filter( !(is.na(SampleID)))

pdf("../../figures/luminex_pca_clad_status.pdf")
ggplot(plot_pca, aes(x = PC1, y = PC2, color = as.factor(status))) + 
    theme_classic() + 
    geom_point() + 
    stat_ellipse() 
ggplot(plot_pca, aes(x = PC3, y = PC4, color = as.factor(status))) + 
    theme_classic() + 
    geom_point() + 
    stat_ellipse() 
dev.off()

summary(aov(PC1 ~ status, data = plot_pca))
summary(aov(PC2 ~ status, data = plot_pca))

## permanova, remove healthy
transplant  = merged
rownames(transplant) = transplant$alt_id
transplant = transplant[ !(is.na(transplant$status)),]

transplant_pca = luminex_pca$x[ rownames(luminex_pca$x) %in% rownames(transplant),]

## eucl distance via pca
luminex_dist = dist(transplant_pca)
transplant = transplant[rownames(as.matrix(luminex_dist)),]

permanova = adonis(luminex_dist~ as.factor(status), data = transplant)
print(permanova$aov.tab)

## ------------- Association between IP-10 and strep-prev ratio -----------
genus = read.csv("../../data/BAL_genus_filtered.csv", row.names = 1)
IP10_filt = IP10[ IP10$SampleID %in% rownames(genus),]
genus = genus[ IP10_filt$SampleID,]
percent = genus / rowSums(genus)
log_genus = log10(genus + 1)

IP10_filt$prev_strep_log_ratio = log_genus[, "g__Streptococcus"] - log_genus[, "g__Prevotella"]
IP10_filt$IP10_log10 = log10(1 + IP10_filt$Concentration)

cor.test( IP10_filt$prev_strep_log_ratio, IP10_filt$IP10_log10, method = "pearson")
cor.test( IP10_filt$prev_strep_log_ratio, IP10_filt$Concentration, method = "pearson")

pdf("../../figures/IP10_prev_strep.pdf", height = 6, width = 6)
ggplot(IP10_filt, aes(x = prev_strep_log_ratio, y = IP10_log10)) + 
    geom_point() + 
    stat_smooth(method = "lm") + 
    theme(text = element_text(size = 16)) + 
    theme_classic() +
    ylab("log10 IP10 Concentration") +
    xlab("Streptococcus / Prevotella log ratio")
dev.off()

summary(coxph(Surv(clad_free_time, status) ~ prev_strep_log_ratio + IP10_log10, data = IP10_filt))
summary(coxph(Surv(clad_free_time, status) ~ prev_strep_log_ratio + Concentration, data = IP10_filt))


