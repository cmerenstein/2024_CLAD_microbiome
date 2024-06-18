library(tidyverse)
library(vegan)

## get the UniFrac distance between each timepont and the last sample from the same subject
## Then get total distance / total time

meta = read.csv("../data/BAL_meta_filtered.csv", row.names = 1)
meta$Time = as.Date(meta$DateCollected, format = "%m/%d/%Y")
weighted = readRDS("../data/weighed_unifrac_no_contam.rds")

## only get the overlaps - the unifract distance will drop some more samples because
## we removed the contaminants and now more are <1000 reads
meta = meta[ meta$SampleID %in% rownames(weighted),]

## ----- Distance from one's own recipient timepoint ---------------------------
transplants = meta[ !(meta$timepoint %in% c("H", "D")),] ## not healthy or donor

w_last = group_by(transplants, SubjectID) %>%
    arrange(Time) %>%
    mutate(last_sample = lag(SampleID), last_time = lag(Time)) %>%
    ungroup() %>%
    as.data.frame()

## remove the first samples b/c they don't have a prior sample
w_last = w_last[ !(is.na(w_last$last_sample)), ]
w_last$time_diff = w_last$Time - w_last$last_time

## Remove timepoints that are at the same time -- these are errors
w_last_filter = filter(w_last, time_diff != 0)

distance = w_last_filter[, c("SampleID", "SubjectID", "Time", "last_sample",
                            "last_time", "time_diff")]
distance$distance = apply(distance, 1, function(row){
                        d = weighted[ row["SampleID"], row["last_sample"]]
                        return(d) }) 
weighted_no_diag = weighted
diag(weighted_no_diag) <- NA
mean(weighted_no_diag, na.rm = T)
mean(distance$distance)
mean(distance$distance) / mean(weighted_no_diag, na.rm = T)
sum(distance$distance > mean(weighted_no_diag, na.rm = T)) / nrow(distance)

pdf("../figures/distance_by_time.pdf")
ggplot(filter(distance, time_diff < 500), aes(x = time_diff, y = distance)) + 
    geom_point() + 
    theme_classic() 
dev.off()

cumulative_distance = group_by(distance, SubjectID) %>%
        arrange(Time) %>%
        mutate( c_dist = cumsum(distance), c_time = cumsum(as.numeric(time_diff)),
                n_timepoints = n()) %>%
        ungroup() %>%
        arrange(SubjectID, Time) %>%
        as.data.frame() 

pdf("../figures/cumulative_distance.pdf")

filter(cumulative_distance, c_time < 600) %>%
        group_by(SubjectID) %>% slice_max(order_by = c_time) %>%
        ungroup() %>%
ggplot(  aes(x = c_time, y = c_dist, color = as.factor(n_timepoints))) + 
    geom_point() +
    stat_smooth( method = "lm") + 
    theme_classic() + 
    xlab("Cumulative Sampling Time") +
    ylab("Cumulative Weighted UniFrac Distance")
dev.off()

### ----------------------- CoxPH model with distance per timepoint --------------
library(tidyverse)
library(survival)
library(survminer)

last_tp = filter(cumulative_distance, c_time < 600) %>%
        group_by(SubjectID) %>% slice_max(order_by = c_time) %>%
        ungroup() %>% 
        mutate(dist_per_sample = c_dist / n_timepoints) %>%
        as.data.frame()

## ---- Read in time to clad & censor data -----
KM = read.csv("../data/merged_KM_data.csv")
KM$microbiome_pid = as.character(KM$microbiome_pid)
KM = KM[,c("microbiome_pid", "clad_free_time", "status")]

merged = inner_join(KM, last_tp, by = c("microbiome_pid" = "SubjectID")) %>% unique()
rownames(merged) = merged$SampleID

merged$dist_quartile = merged$dist_per_sample %>% cut(3, labels = c("low", "mid", "high"))

fit1 = survfit(Surv(clad_free_time, status)~ dist_quartile, data = merged)

#pdf("../figures/survival_distance_traveled.pdf")
#ggsurvplot(fit1, conf.int=TRUE, pval=TRUE)
#dev.off()


coxph(Surv(clad_free_time, status) ~ dist_per_sample, data = merged)








