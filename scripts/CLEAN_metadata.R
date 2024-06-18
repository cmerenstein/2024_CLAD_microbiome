library(tidyverse)

## %%%%%%%%%%%%%%%%%%%    Load and Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## ------------- Read in sample metadata as was given to me by the clinical team ----
meta = read.csv("../data/old_metadata.csv")
meta$timepoint = ifelse( grepl("HUP", meta$SampleID), "H", meta$timepoint) # healthy controls 
rownames(meta) = meta$original.sample.ids
meta$better_id = paste(meta$SubjectID, meta$timepoint, meta$SampleType, sep = ".")
print(dim(meta))


## --------------- Fix a few known errors in the metadata ----------------
meta = filter(meta, diagnosis != "CF")
meta = filter(meta, SubjectID != "641001019") ## actually a CF patient, shouldn't be here
meta$diagnosis[meta$diagnosis == "Non-IPF ILD"] <- "Non-IPF"
meta$diagnosis[meta$diagnosis == "OB pre-transplant"] <- "BOS pre-txp"
meta$diagnosis = ifelse(meta$SubjectID == "641001007", "COPD", meta$diagnosis)
meta$transplant = ifelse(meta$SubjectID == "641000911", "BLT", meta$transplant)
meta$transplant = ifelse(meta$SubjectID == "641000072", "UNK", meta$transplant)
meta$timepoint[meta$timepoint == "RT"] <- "R"

meta$transplant_date = gsub("2006", "2016", meta$transplant_date)

## --------- Load genus data to collapse ------------
genus = read.csv("../qiime_data/ASV/counts_taxa_level_5.csv", row.names = 1)
genus = genus[rownames(genus) %in% rownames(meta),]
genus_df = data.frame(sample_rep = rownames(genus), genus)
print(dim(genus_df))

## ------------------ combine replicates ----------------------------
genus_df = data.frame(sample = rownames(genus), genus)
tidy = gather(genus_df, "genus", "count", -sample)
tidy$ID = meta[tidy$sample, "SampleID"]

tidy_collapsed = group_by(tidy, ID, genus) %>% summarize(count = sum(count)) %>%
    ungroup() %>% as.data.frame()

collapsed = spread(tidy_collapsed, "genus", "count", fill = 0)
collapsed_mat = collapsed[, 2:ncol(collapsed)]
rownames(collapsed_mat) = collapsed$ID

print(dim(collapsed_mat))

## ------  collapse metadata by replicate ----------------
meta_collapsed = meta[, colnames(meta) != "original.sample.ids"] %>% unique()
rownames(meta_collapsed) = meta_collapsed$SampleID
print(dim(meta_collapsed))

## ---- Get just the BAL samples ----------------
BAL = meta_collapsed[ meta_collapsed$SampleType == "BAL fluid",]
BAL_depth = rowSums(collapsed_mat)[rownames(BAL)]
BAL$depth = BAL_depth

BAL = BAL[ BAL$SampleID %in% rownames(collapsed_mat),]
rownames(BAL) = BAL$SampleID
print(dim(BAL))

## ----- Handle rerun samples ---------------------------
## For now, going to remove run 4 samples that have replicates in run 6
BAL$better_id = paste(BAL$SubjectID, BAL$timepoint, BAL$SampleType, sep = '.')

duplicates = group_by(BAL, better_id) %>% summarize(n = n() ) %>% ungroup() %>% data.frame()
duplicates = duplicates[duplicates$n == 2, "better_id"]

BAL = filter(BAL, !(better_id %in% duplicates & run == 4))
print(dim(BAL))

## --- Remove re-run samples from the full dataframe too (i.e. the one with controls) ---
meta_unfiltered_w_controls = filter(meta_collapsed, SampleType != "BAL fluid" | 
                                        SampleID %in% BAL$SampleID ) %>%
                                filter(SampleID %in% rownames(collapsed_mat))
genus_unfiltered_w_controls = collapsed_mat[ meta_unfiltered_w_controls$SampleID,]

write.csv(genus_unfiltered_w_controls, "../data/genus_unfiltered_w_controls.csv")
write.csv(meta_unfiltered_w_controls, "../data/meta_unfiltered_w_controls.csv")

## --- Write unfiltered and filtered genus dataframe --------
taxa_BAL = collapsed_mat[ BAL$SampleID,]
print(dim(taxa_BAL))
write.csv(taxa_BAL, "../data/BAL_genus_unfiltered.csv")
write.csv(BAL, "../data/BAL_meta_unfiltered.csv")

threshold_row = 1000
filtered_row = taxa_BAL[rowSums(taxa_BAL) >= threshold_row,]

threshold_col = 100 # at least 100 reads in any one sample
filtered = filtered_row[, colSums(filtered_row > threshold_col) >= 1]
write.csv(filtered, "../data/BAL_genus_filtered.csv")

BAL_filtered = BAL[rownames(filtered),] 
write.csv(BAL_filtered,"../data/BAL_meta_filtered.csv")


## %%%%%%%%%%%% Load CLAD data and get time to CLAD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLAD = read.csv("../data/2024-04-19microbiomecladdata.csv")
ltog_to_mb = CLAD[, c("ltog_id", "microbiome_pid")]

### 
CLAD$microbiome_pid  = as.character(CLAD$microbiome_pid)
CLAD = CLAD[ CLAD$microbiome_pid %in% c(BAL_filtered$SubjectID),] 

CLAD$pft_date = as.Date(CLAD$pft_date, format = "%m/%d/%Y")
CLAD$transplant_date = as.Date(CLAD$transplant_date, format = "%m/%d/%Y")

end_CLAD = group_by(CLAD, microbiome_pid) %>% filter(pft_date == max(pft_date)) %>% 
        ungroup() %>% as.data.frame() %>% filter(type_clad == "Definite")

end_CLAD$CLAD_onset = end_CLAD$pft_date - end_CLAD$clad_duration
end_CLAD$time_to_CLAD = as.numeric(end_CLAD$CLAD_onset - end_CLAD$transplant_date)

time_to_clad = end_CLAD[, c("microbiome_pid", "CLAD_onset", "time_to_CLAD")]

merged_CLAD = merge(CLAD, time_to_clad, by = "microbiome_pid", all.x = T)
clean_CLAD = merged_CLAD[, c("microbiome_pid", "transplant_id", "transplant_date",
                            "pft_date", "days_from_transplant", "CLAD_onset", "time_to_CLAD", 
                            "clad_duration", "type_clad")]
CLAD_by_patient = group_by(clean_CLAD, microbiome_pid) %>%
                    top_n(1, days_from_transplant) %>%
                    ungroup() %>% as.data.frame()
CLAD_by_patient$clad_free_time = ifelse(is.na(CLAD_by_patient$time_to_CLAD), 
                                    CLAD_by_patient$days_from_transplant, CLAD_by_patient$time_to_CLAD)
CLAD_by_patient$status = as.numeric( !(is.na(CLAD_by_patient$time_to_CLAD)))

## before writing, add in relevant demographic information
demo = read.csv("../data/2024-05-03additionalmicrodata.csv")
demo = demo[,c("microbiome_pid", "age", "gender", "BMI")]
demo_CLAD = merge(CLAD_by_patient, demo, by = "microbiome_pid")

write.csv(demo_CLAD, "../data/merged_KM_data.csv")


## %%%%%%%%%% Get the phylum level data cleaned up too %%%%%%%%%%%%%%%%%%%%%%%%%%%
phylum = read.csv("../qiime_data/ASV/counts_taxa_level_1.csv", row.names = 1)
phylum = phylum[rownames(phylum) %in% rownames(meta),]
phylum_df = data.frame(sample_rep = rownames(phylum), phylum)
print(dim(phylum_df))

## ------------------ combine replicates ----------------------------
phylum_df = data.frame(sample = rownames(phylum), phylum)
tidy = gather(phylum_df, "phylum", "count", -sample)
tidy$ID = meta[tidy$sample, "SampleID"]

tidy_collapsed = group_by(tidy, ID, phylum) %>% summarize(count = sum(count)) %>%
    ungroup() %>% as.data.frame()

collapsed = spread(tidy_collapsed, "phylum", "count", fill = 0)
collapsed_mat = collapsed[, 2:ncol(collapsed)]
rownames(collapsed_mat) = collapsed$ID

print(dim(collapsed_mat))
phylum_BAL = collapsed_mat[ BAL_filtered$SampleID,]
write.csv(phylum_BAL, "../data/BAL_collapsed_phylum_filtered.csv")





