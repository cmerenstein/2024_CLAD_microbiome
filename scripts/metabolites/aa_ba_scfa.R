library(tidyverse)
library(vegan)
library(ape)
library(survival)
library(survminer)

meta = read.csv("../../data/BAL_meta_unfiltered.csv", row.names = 1)

KM = read.csv("../../data/merged_KM_data.csv")[, c("microbiome_pid", "status", "age", "clad_free_time")]
KM$microbiome_pid = as.character(KM$microbiome_pid)
rownames(KM)  = KM$microbiome_pid

metabolomics = read.csv("../../data/metabolites/aa_ba_scfa_metabolomics_run1.csv")
metabolomics$SubjectID = as.character(metabolomics$SubjectID)

## seems to be one duplicate subject ID. One doesn't correspond to any proper sample and is a mistake
metabolomics = filter(metabolomics, !(SubjectID =="641000838" & CollectionDate == "2/6/2017"))
rownames(metabolomics) = metabolomics$SubjectID

## samples_with_info
intersect_samples = intersect(rownames(KM), rownames(metabolomics))
KM = KM[intersect_samples,]
metabolomics = metabolomics[intersect_samples, ]

## ----- get in matrix form -----------
metabolomics_mat = as.matrix(metabolomics[,3:ncol(metabolomics)])

## NA as zeros?
metabolomics_mat[is.na(metabolomics_mat)] <- 0

dist = vegdist(metabolomics_mat, "bray")
dist_mat = as.matrix(dist)[rownames(KM), rownames(KM)]
adonis(dist_mat~status, data = KM)$aov.tab

## inidividual mets
metabolomics_filter = metabolomics_mat[, colSums(metabolomics_mat > 0) > 10]

met_coxPH = lapply(colnames(metabolomics_filter), function(met){
                   KM_copy = KM
                   KM_copy$met = metabolomics_filter[rownames(KM_copy), met] 
                   cox = summary(coxph(Surv(clad_free_time, status)~met + age, data = KM_copy))
             return(data.frame( metabolte = met, HR = cox$coefficients[1, 2], p = cox$coefficients[1,5]))
            })
met_coxPH_df = do.call("rbind", met_coxPH)
met_coxPH_df$FDR = p.adjust(met_coxPH_df$p, "BH")
arrange(met_coxPH_df, p)
write.csv(arrange(met_coxPH_df, p), "aa_scfa_coxph.csv", row.names = F, quote = F)



