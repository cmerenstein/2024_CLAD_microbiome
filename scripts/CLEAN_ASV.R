library(tidyverse)
library(phyloseq)
library(ape)

tree <- read.tree("../qiime_data/tree.nwk")
tree <- ape::multi2di(tree) ## need to fix the tree after pruning contaminats


## use the old metadata to be able to collapse samples
meta = read.csv("../data/old_metadata.csv")
meta$timepoint = ifelse( grepl("HUP", meta$SampleID), "H", meta$timepoint)
rownames(meta) = meta$original.sample.ids
meta$better_id = paste(meta$SubjectID, meta$timepoint, meta$SampleType, sep = ".")
print(dim(meta))

## ------ Load ASV data and collapse replicates -----------
asv <- read.table("../qiime_data/DADA2/feature_table_clean.txt",
                    header = TRUE, row.names = 1, sep = "\t") %>% t()
asv = asv[rownames(asv) != "Undetermined",]

## filter super rare ASVs to make everything take less time
asv = asv[, colSums(asv) > 10]

asv_df = data.frame(sample = rownames(asv), asv)
colnames(asv_df) = gsub("X", "", colnames(asv_df))
tidy = gather(asv_df, "asv", "count", -sample)
tidy$ID = meta[tidy$sample, "SampleID"]
tidy = filter(tidy, !(is.na(ID)))

tidy_collapsed = group_by(tidy, ID, asv) %>% summarize(count = sum(count)) %>%
    ungroup() %>% as.data.frame()

collapsed = spread(tidy_collapsed, "asv", "count", fill = 0)
collapsed_mat = collapsed[, 2:ncol(collapsed)]
rownames(collapsed_mat) = collapsed$ID

print(dim(collapsed_mat))

write.csv(collapsed_mat, "../data/ASV_unfiltered.csv")

## ------------- get just BAL samples and handle re-runs. 
BAL = read.csv("../data/BAL_meta_unfiltered.csv", row.names = 1)
asv_BAL = collapsed[ collapsed$ID %in% rownames(BAL),] 
write.csv(asv_BAL, "../data/BAL_asv_unfiltered.csv", row.names = F)

read_asv = read.table("../data/BAL_asv_unfiltered.csv", sep = ',', header = T, check.names = F)

