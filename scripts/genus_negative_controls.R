library(tidyverse)
library(vegan)
library(ape)

## ---------- Read in all metadata 
meta = read.csv("../data/old_metadata.csv")
rownames(meta) = meta$original.sample.ids

genus = read.csv("../qiime_data/ASV/counts_taxa_level_5.csv", row.names = 1)
genus = genus[rownames(genus) %in% rownames(meta),]

## ----- Filter very loosely @ > 100 counts -------------
genus = genus[rowSums(genus) > 100,]
meta = meta[rownames(genus),]

percent = genus / rowSums(genus)

### ------------- Plot genera abundance in controls  -------------------------
controls = meta[ meta$SampleType != "BAL fluid", ]
controls$run = as.factor(controls$run)

genus_controls = genus[ rownames(controls), ]
genus_controls = genus_controls[ , names(tail(sort(colSums(genus_controls)), n = 40))]

genus_control_taxa = genus[, colnames(genus_controls)]

genus_df = data.frame(sample = rownames(genus_control_taxa), genus_control_taxa)
tidy = gather(genus_df, "taxa", "count", -sample)
tidy$type = meta[ tidy$sample, "SampleType"]
tidy$control = ifelse(tidy$type != "BAL fluid", "Control", "BAL")
tidy$batch = as.factor( meta[ tidy$sample, "run"])

pdf("../figures/contamination_genera.pdf", width = 10, height = 20)

ggplot( tidy, aes( x = batch, y = count, fill = control) ) + 
    geom_boxplot() + 
    facet_wrap(~taxa, scales = "free", ncol = 4) + 
    theme_classic()
dev.off()

## --------------- genera that are more abundant in negative controls ---------
## use collapsed samples and metadata
meta_collapsed = read.csv("../data/meta_unfiltered_w_controls.csv", row.names = 1)
genus_collapsed = read.csv("../data/genus_unfiltered_w_controls.csv", row.names = 1)

tidy_no_contam_list = lapply( unique(meta_collapsed$run), function(run_n){

    this_run = filter(meta_collapsed, run == run_n)
    genus_this_run = genus_collapsed[ rownames(this_run),]
    
    genus_df = data.frame(sample = rownames(genus_this_run), genus_this_run)
    tidy_this_run = gather(genus_df, "taxa", "value", -sample)
    tidy_this_run$BAL = this_run[ tidy_this_run$sample, "SampleType" ] == "BAL fluid"
   
    is_contam_list = lapply( unique(tidy_this_run$taxa), function(tax){

        tidy_tax = filter(tidy_this_run, taxa == tax)
        pval = wilcox.test( value ~ BAL, data = tidy_tax)$p.value
        diff = mean( tidy_tax[ tidy_tax$BAL, "value"]) -
                    mean(tidy_tax[ !(tidy_tax$BAL), "value"]) ## sample - control
        return( data.frame( run = run_n, taxa = tax, contam = (pval < 0.05 & diff < 0)) )
    }) 

    is_contam = do.call("rbind", is_contam_list)
    print(filter(is_contam, contam))
    contam_df = full_join(tidy_this_run, is_contam, by = "taxa")
    
    contam_df$value = ifelse( contam_df$contam, 0, contam_df$value)
    
    return(contam_df)
})
tidy_no_contam = do.call("rbind", tidy_no_contam_list)
genus_no_contam = tidy_no_contam[, c("sample", "taxa", "value")] %>%
                    spread(taxa, value)
no_contam_mat = genus_no_contam[,2:ncol(genus_no_contam)]
rownames(no_contam_mat) = genus_no_contam$sample

write.csv(genus_no_contam, "../data/genus_no_contam.csv", row.names = F)



