library(tidyverse)
library(vegan)

## ----------------- Load, filter, combine rare taxa -------------
meta = read.csv("../data/BAL_meta_filtered.csv", row.names = 1)
meta$timepoint[ meta$timepoint == "RT" ] <- "R"
genus = read.csv( "../data/genus_no_contam.csv", row.names = 1)
genus = genus[rownames(meta), colSums(genus > 10) > 0]
dim(genus)

## combine rare and common taxa
rare_taxa = genus[, colSums( genus > 10 ) < (0.1 * ncol(genus))]
combined_rare_taxa = rowSums(rare_taxa)
common_taxa = genus[, colSums(genus > 10) >= (0.1 * ncol(genus))]

genus_reduced = cbind(common_taxa, combined_rare_taxa)

##  now have to re-filter because some no longer rowSums to 1000
genus_reduced = genus_reduced[ rowSums(genus_reduced) >= 1000,]
meta = meta[rownames(genus_reduced),]
percent = genus_reduced / rowSums(genus_reduced)


## ---------------- Genus heatmap -------------------------------
library(pheatmap)

get_order <- function(samples) { return(hclust(dist(percent[ rownames(percent) %in% samples,]))$order)}
get_order(meta[ meta$timepoint == "R", "SampleID"])

sample_order = lapply(c("H", "D", "R", "DC", "6W", "3M", "6M", "12M"), function(tp){
    samples = meta[meta$timepoint == tp, "SampleID"]
    reordered = samples[get_order(samples)]
    return(data.frame(SampleID = reordered, timepoint = tp))
})

shpitul = do.call("rbind", sample_order[1:4])
rownames(shpitul) = shpitul$SampleID
post = do.call("rbind", sample_order[5:8])
rownames(post) = post$SampleID

percent_shpitul = percent[shpitul$SampleID,]
percent_post = percent[post$SampleID, ]

## filter for plotting
percent_shpitul = percent_shpitul[ , colSums( percent_shpitul > .3) > 0]
percent_post = percent_post[ , colSums( percent_post > .3) > 0]

pdf("../figures/no_contam_genus_heatmap.pdf", width = 12, height = 8)

pheatmap( t(percent_shpitul), annotation_col = shpitul[,c("timepoint"), drop = F],
            labels_col = character(nrow(percent_shpitul)), cluster_cols = F,
            treeheight_row = 0 )
pheatmap( t(percent_post), annotation_col = post[,c("timepoint"), drop = F],
            labels_col = character(nrow(percent_post)), cluster_cols = F,
            treeheight_row = 0 )
dev.off()

