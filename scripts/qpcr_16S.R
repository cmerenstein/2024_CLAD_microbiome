library(tidyverse)
library(survival)
library(survminer)

## ------------------------- read the metadata and merge -------------------------------
KM = read.csv("../data/merged_KM_data.csv", row.names = 1)
meta = read.csv("../data/BAL_meta_filtered.csv", row.names = 1)
meta$timepoint[ meta$timepoint == "RT"] <- "R"
transplant = meta[ !(meta$timepoint %in% c("H", "D")),]

merged = merge(KM, transplant, by.x = "microbiome_pid", by.y = "SubjectID")
sample_id_date = transplant[,c("SampleID", "SubjectID", "DateCollected", "Description")]

## -------------------------- read qpcr and rerun qpcr --------------------------------

reruns = read.csv("../data/qPCR/reruns_copy_number.csv")
copies = read.csv("../data/qPCR/Copies_ID.csv") %>%
    merge(reruns, by = "SampleID", all = T)
copies = copies[!(is.na(copies$Copy.)),]

copies$copies = ifelse(is.na(copies$Copy..Re.Run), copies$Copy., copies$Copy..Re.Run)

copies = copies[,c("SampleID", "copies")]
copies$log10_Copies = log10(1 + copies$copies)


## --------------------------- magic happens in this jawn ------------------
qpcr = merge( copies, transplant, by = "SampleID")

coxph = lapply( unique(qpcr$timepoint), function(tp){
    qpcr_filt = qpcr[ qpcr$timepoint == tp,] %>%
                merge(KM, by.x = "SubjectID", by.y = "microbiome_pid")
    cox = coxph(Surv(clad_free_time, status) ~ log10_Copies + age, data = qpcr_filt)
    return(data.frame(tp, p = summary(cox)$coefficients[1,5],
                        OR = summary(cox)$coefficient[1,2])) })
coxph = do.call("rbind", coxph)
coxph$FDR = p.adjust(coxph$p, "BH")

## ----- merge and plot --------------
to_plot = merge(KM, qpcr, by.x = "microbiome_pid", by.y = "SubjectID")
to_plot$timepoint = factor(to_plot$timepoint, c( "R", "DC", "6W", "3M", "6M", "12M"))

pdf("../figures/qpcr_by_CLAD.pdf")
ggplot(to_plot, aes(x = timepoint, fill = as.factor(status), y = log10_Copies)) + 
    geom_boxplot() + 
    theme_classic() + 
    guides(fill = guide_legend(title = "CLAD")) + 
    theme(text = element_text(size = 20))
dev.off()

pvals_wilcox = sapply(unique(qpcr$timepoint), function(tp){
                filt = to_plot[to_plot$timepoint == tp,]
                wilcox = wilcox.test( log10_Copies ~ status, data = filt)
                return(wilcox$p.value)
        })
p.adjust(pvals_wilcox, "BH")

## ------- plot with continuous time -----------
qpcr$time_from_transplant = as.numeric( as.Date(qpcr$DateCollected, "%m/%d/%Y") -
                            as.Date(qpcr$transplant_date, "%m/%d/%Y"))
qpcr_w_time = merge(KM, qpcr, by.x = "microbiome_pid", by.y = "SubjectID") %>%
                filter(time_from_transplant < 400 & time_from_transplant >= 0)
qpcr_w_time$CLAD_status = c("Non-CLAD", "CLAD")[qpcr_w_time$status + 1]

pdf("../figures/qpcr_time_from_transplant.pdf")
ggplot(qpcr_w_time, aes(x = time_from_transplant, y = log10_Copies, color = CLAD_status)) +
    theme_classic() +
    geom_point() +
#    geom_line(aes(group = microbiome_pid), alpha = 0.2, linewidth = .5) +
    scale_color_manual(values = c("#f75965", "#a2d6fa")) +
    stat_smooth() +
    xlab("Time from transplant") +
    ylab("16S gene log10 copies per mL") +
    guides(color = guide_legend(title = "CLAD status")) +
    theme(text = element_text(size = 16))
dev.off()







