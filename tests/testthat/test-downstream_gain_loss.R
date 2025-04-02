library(testthat)
library(phylosuite)
library(dplyr)
library(stringr)
library(ggplot2)
library(ape)
tr <- phylosuite::tr
df <- phylosuite::df
geno <- phylosuite::genotype_mat
df <- left_join(df,geno %>% mutate(tip_name_var = rownames(.)))
df$isolate_no <- df$tip_name_var
pheno <-  'colistin_ns'
geno <-  colnames(geno)[1]
tip_name_var <- "tip_name_var"

asr_pheno_obj <- asr(df = df,tr = tr,tip_name_var = tip_name_var ,pheno = pheno,model="ER",node_states = "joint",conf_threshold=NULL)
asr_geno_obj <- asr(df = df,tr = tr,tip_name_var = tip_name_var ,pheno = pheno,model="ER",node_states = "joint",conf_threshold=NULL)
asr_downstream <- downstream_gain_loss(comparitor_parent_child_df = asr_geno_obj$parent_child_df,trait_parent_child_df = asr_pheno_obj$parent_child_df,tr = tr,node_states = "joint",confidence = NULL)

expect_s3_class(asr_downstream, "data.frame")

expect_equal(object = nrow(asr_downstream),expected = 1)
expect_equal(object = ncol(asr_downstream),expected = 14)
