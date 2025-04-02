library(testthat)
library(phyloAMR)
library(dplyr)
library(stringr)
library(ggplot2)
library(ape)
tr <- phyloAMR::tr
df <- phyloAMR::df
pheno <-  'colistin_ns'
tip_name_var <- "tip_name_var"

asr_out <- phyloAMR::asr(df = df,tr = tr,tip_name_var = tip_name_var,pheno = pheno,model="ER",node_states="joint")

expect_s3_class(asr_out$corHMM_out, "corhmm")
expect_s3_class(asr_out$parent_child_df, "data.frame")

expect_equal(object = nrow(asr_out$parent_child_df),expected = nrow(tr$edge))

expect_error(phyloAMR::asr(df,tr,tip_name_var,pheno,model="ER",node_states="t",conf_threshold=0.875))

expect_error(phyloAMR::asr(df,tr,tip_name_var,pheno,model="RDS",node_states="t",conf_threshold=0.875))
