# Load necessary libraries
library(ape)
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)

# Replace 'path_to_tree' with the path to your newick-formatted tree file
tree <- read.tree("../qiime_data/tree.nwk")
tree <- ape::multi2di(tree) ## need to fix the tree after pruning contaminats

## Just the BAL filtered samples
BAL = read.csv("../data/BAL_meta_filtered.csv", row.names = 1)

asv = readRDS("../data/BAL_asv_unfiltered.rds")
rownames(asv) = asv$ID
asv = asv[rownames(BAL), ]
asv_mat = asv[,2:ncol(asv)]
asv_mat = asv_mat[,colSums(asv_mat) > 100]

## rarefied / just normalize to cpm
rare = asv_mat / rowSums(asv_mat) * 100000
rare_phylo = phyloseq( otu_table(rare, taxa_are_rows = F), tree)
rare_weighted = UniFrac(rare_phylo, weighted = T)


## ---------- PCoA weighted UniFrac -----------------------------------
## load genus data for fiting
set.seed(19143)
ord = pcoa(as.matrix(rare_weighted))
vectors = ord$vectors

pcoa_df = data.frame(sample = rownames(vectors), Axis.1 = vectors[,1],
                     Axis.2 = vectors[,2], Axis.3 = vectors[,3], Axis.4 = vectors[,4],
                     timepoint = BAL$timepoint,
                     run = as.factor(BAL$run)) %>%
                        filter( !(is.na(timepoint)))

pdf("../figures/pcoa_batch_effect_WU.pdf")
ggplot( pcoa_df, aes( x = Axis.1, y = Axis.2, color = run)) +
        theme_classic() +
        geom_point() +
        stat_ellipse() +
        ggtitle("pcoa")

ggplot( pcoa_df, aes( x = Axis.3, y = Axis.4, color = run)) +
        theme_classic() +
        geom_point() +
        stat_ellipse() +
        ggtitle("pcoa")

dev.off()
 



