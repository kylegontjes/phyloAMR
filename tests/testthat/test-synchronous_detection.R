library(testthat)
library(phyloAMR)
library(dplyr)
library(stringr)
library(ggplot2)
library(ape)
tr <- phyloAMR::tr
df <- phyloAMR::df
geno <- phyloAMR::genotype_mat
df <- left_join(df,geno %>% mutate(tip_name_var = rownames(.)))
df$isolate_no <- df$tip_name_var
pheno <-  'colistin_ns'
geno <-  colnames(geno)[1]
tip_name_var <- "tip_name_var"

asr_pheno_obj <- asr(df = df,tr = tr,tip_name_var = tip_name_var ,pheno = pheno,model="ER",node_states = "joint",conf_threshold=NULL)
asr_geno_obj <- asr(df = df,tr = tr,tip_name_var = tip_name_var ,pheno = pheno,model="ER",node_states = "joint",conf_threshold=NULL)
asr_sync <- synchronous_detection(comparitor_parent_child_df = asr_geno_obj$parent_child_df,trait_parent_child_df = asr_pheno_obj$parent_child_df)

expect_s3_class(asr_sync, "data.frame")

expect_equal(object = nrow(asr_sync),expected = 1)
expect_equal(object = ncol(asr_sync),expected = 14)
