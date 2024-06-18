library(tidyverse)
library(vegan)
library(ape)
library(phyloseq)


## ------ Is an asv more abundant in controls vs samples? ---------
## ------ On a per-batch basis
## ------ Use collapsed samples and metadata
meta_collapsed = read.csv("../data/meta_unfiltered_w_controls.csv", row.names = 1)
#asv_collapsed = read.csv("../data/ASV_unfiltered.csv", row.names = 1)
asv_collapsed = readRDS("../data/ASV_unfiltered.rds")
asv_collapsed = asv_collapsed[ rownames(meta_collapsed),]

tidy_no_contam_list = lapply( unique(meta_collapsed$run), function(run_n){

    this_run = filter(meta_collapsed, run == run_n)
    asv_this_run = asv_collapsed[ rownames(this_run),]
    asv_this_run = asv_this_run[, colSums(asv_this_run) > 100,]
    
    asv_df = data.frame(sample = rownames(asv_this_run), asv_this_run)
    tidy_this_run = gather(asv_df, "taxa", "value", -sample)
    tidy_this_run$BAL = this_run[ tidy_this_run$sample, "SampleType" ] == "BAL fluid"
   
    is_contam_list = lapply( unique(tidy_this_run$taxa), function(tax){

        tidy_tax = filter(tidy_this_run, taxa == tax)
        pval = wilcox.test( value ~ BAL, data = tidy_tax)$p.value
        diff = mean( tidy_tax[ tidy_tax$BAL, "value"]) -
                    mean(tidy_tax[ !(tidy_tax$BAL), "value"]) ## sample - control
        return( data.frame( run = run_n, taxa = tax, contam = (pval < 0.05 & diff < 0)) )
    }) 

    is_contam = do.call("rbind", is_contam_list)
    contam_df = full_join(tidy_this_run, is_contam, by = "taxa")
    
    contam_df$value = ifelse( contam_df$contam, 0, contam_df$value)
    
    return(contam_df)
})
tidy_no_contam = do.call("rbind", tidy_no_contam_list)
asv_no_contam = tidy_no_contam[, c("sample", "taxa", "value")] %>%
                    spread(taxa, value, fill = 0)
no_contam_mat = asv_no_contam[,2:ncol(asv_no_contam)]
rownames(no_contam_mat) = asv_no_contam$sample

saveRDS(asv_no_contam, "../data/asv_no_contam.rds")
#write.csv(asv_no_contam, "../data/asv_no_contam.csv", row.names = F)

no_contam_mat = asv_no_contam[,2:ncol(asv_no_contam)]
rownames(no_contam_mat) = asv_no_contam$sample

## ----- see how that changes batch effecs ------------
no_contam_filter = no_contam_mat[ rowSums(no_contam_mat) > 1000 & 
                            meta_collapsed[rownames(no_contam_mat), "SampleType"] == "BAL fluid",]

set.seed(19143)
rare = no_contam_filter / rowSums(no_contam_filter)  ## not actually rarefy, just convert to %

tree <- read.tree("../qiime_data/tree.nwk")
tree <- ape::multi2di(tree) ## need to fix the tree after pruning contaminats

colnames(rare) = gsub("X", "", colnames(rare))
rare_phylo = phyloseq( otu_table(rare, taxa_are_rows = F), tree)
rare_weighted = UniFrac(rare_phylo, weighted = T)

saveRDS(as.matrix(rare_weighted), "../data/weighed_unifrac_no_contam.rds")

# -------------------- PCOA plot of UniFrac  by batch --------------
set.seed(19002)
ordinates = pcoa(rare_weighted)

str(ordinates) ## look at relative correlation for variance explained
vectors = ordinates$vectors

pcoa_df = data.frame(sample = rownames(vectors), Axis.1 = vectors[,1],
                    Axis.2 = vectors[,2], Axis.3 = vectors[,3], Axis.4 = vectors[,4],
                    run = as.factor(meta_collapsed[rownames(vectors), "run"]))


pdf("../figures/no_contam_WU_1000.pdf")
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











