library(tidyverse)

KM = read.csv("../data/merged_KM_data.csv")

## format metadata
meta_filt = read.csv("../data/BAL_meta_filtered.csv", row.names = 1)
meta = read.csv("../data/BAL_meta_unfiltered.csv", row.names = 1)
meta$high_quality = rownames(meta) %in% rownames(meta_filt)

meta$time_post_transplant = as.Date(meta$DateCollected, "%m/%d/%Y") - 
                            as.Date(meta$transplant_date, "%m/%d/%Y")
meta = filter(meta, time_post_transplant >= -1)

## R and D are both day 0 but plot them in the negatives to see better
meta$time_to_plot = ifelse(meta$timepoint == "R", -40, 
                        ifelse(meta$timepoint == "D", -80, meta$time_post_transplant))

## load live/dead data
demo = read.csv("../data/2024-05-03additionalmicrodata.csv")
alive = demo[,c("PX_STAT", "PX_STAT_DATE", "microbiome_pid")]
rownames(alive) = as.character(alive$microbiome_pid)

meta$live_dead = alive[meta$SubjectID, "PX_STAT"]
meta$live_dead_date = as.Date(alive[meta$SubjectID, "PX_STAT_DATE"], format = "%d-%b-%y")
meta$live_dead_followup = meta$live_dead_date - as.Date(meta$transplant_date, "%m/%d/%Y")

## -------- Merge and plot -----------------
CLAD = KM[, c("microbiome_pid", "clad_free_time", "days_from_transplant", "status")] %>%
                rename(total_followup = days_from_transplant)
merged = merge(meta, CLAD, by.x = "SubjectID", by.y = "microbiome_pid") %>%
                arrange(status, total_followup)
merged$status = as.factor(merged$status)

## difference in follow up for PFT and status
merged$followup_diff = merged$live_dead_followup - merged$total_followup
merged$followup_status = ifelse(merged$status ==1, "CLAD",
                           ifelse(merged$followup_diff > 90, "non-CLAD",
                            ifelse(merged$live_dead == "D", "non-CLAD death", "non-CLAD")))

filtered = merged %>% filter(time_to_plot < 500) %>% 
            filter(total_followup < 5000) ## can't have this much followup
filtered = arrange(filtered, -clad_free_time)
filtered$SubjectID = factor(filtered$SubjectID, levels = unique(filtered$SubjectID))

pdf("../figures/figure_1_sampling.pdf")
ggplot(filtered, aes(x = time_to_plot, y = SubjectID, shape = high_quality)) + 
    theme_classic() + 
    geom_point() + 
    geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 550) + 
    scale_shape_manual(values = c(4, 1)) + 
    xlab("time post transplant") + 
    theme(axis.text = element_blank()) + 
    theme(axis.text = element_blank()) +
    theme(axis.ticks = element_blank()) +
    xlim(c(-80, 1050)) + 
    geom_segment(x = 550, aes(xend = (clad_free_time / 10) + 550, 
                 y = SubjectID, yend = SubjectID, color = followup_status)) + 
    scale_color_manual(values = c("#f75965", "#a2d6fa", "grey" )) + 
    theme(legend.position = c(.9,.25)) +
    geom_vline(xintercept = (550 + 5 * (365 / 13)), alpha = 0.25, linetype = "dashed") 
    
dev.off()


