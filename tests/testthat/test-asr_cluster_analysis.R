library(testthat)
library(phyloAMR)
library(dplyr)
library(stringr)
library(ggplot2)
library(ape)
tr <- phyloAMR::tr
df <- phyloAMR::df
df$isolate_no <- df$tip_name_variable
trait <-  'colistin_ns'
tip_name_variable <- "tip_name_variable"

asr_obj <- asr(df = df,tr = tr,tip_name_variable = tip_name_variable ,pheno = pheno,model="ARD",node_states = "joint")
asr_cluster <- asr_cluster_detection(df = df,tr = tr,trait = trait,tip_name_variable = "tip_name_variable",patient_id = "PatientID",node_states = 'joint',parent_child_df = asr_obj$parent_child_df, simplify_faux_clusters = FALSE, simplify_revertant = TRUE, collapse_cluster = TRUE)
asr_cluster_analysis <- asr_cluster_analysis(tip_data_df = asr_cluster)

expect_s3_class(asr_cluster_analysis, "data.frame")

expect_equal(object = nrow(asr_cluster_analysis),expected = 1)
expect_equal(object = ncol(asr_cluster_analysis),expected = 17)
