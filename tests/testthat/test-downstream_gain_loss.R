library(testthat)
library(phyloAMR)
library(dplyr)
library(stringr)
library(ggplot2)
library(ape)
tr <- phyloAMR::tr
df <- phyloAMR::df
geno <- phyloAMR::genotype_mat
geno$tip_name_var <- rownames(geno)
df <- left_join(df,geno)
trait <-  'colistin_ns'
geno <-  colnames(geno)[1]
tip_name_variable <- "tip_name_var"

asr_trait_obj <- asr(df = df,tr = tr,tip_name_variable = tip_name_variable ,trait = trait,model="ER",node_states = "joint", confidence_threshold = NULL)
asr_geno_obj <- asr(df = df,tr = tr,tip_name_variable = tip_name_variable ,trait = trait,model="ER",node_states = "joint", confidence_threshold = NULL)
asr_downstream <- downstream_transitions(comparitor_parent_child_df = asr_geno_obj$parent_child_df, trait_parent_child_df = asr_trait_obj$parent_child_df, tr = tr, node_states = "joint", confidence = NULL)

expect_s3_class(asr_downstream, "data.frame")

expect_equal(object = nrow(asr_downstream),expected = 1)
expect_equal(object = ncol(asr_downstream),expected = 19)
