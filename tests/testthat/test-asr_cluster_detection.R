library(testthat)
library(phyloAMR)
library(dplyr)
library(stringr)
library(ggplot2)
library(ape)
tr <- phyloAMR::tr
df <- phyloAMR::df
trait <-  'colistin_ns'
tip_name_variable <- "tip_name_var"

asr_obj <- asr(df = df,tr = tr,tip_name_variable = tip_name_variable, trait = trait, model="ER", node_states = "joint")
asr_cluster <- asr_cluster_detection(df = df, tr = tr, tip_name_variable = tip_name_variable, patient_id = "PatientID", node_states = 'joint', parent_child_df = asr_obj$parent_child_df, confidence = NULL,simplify_faux_clusters = FALSE, simplify_revertant = TRUE, collapse_cluster = TRUE)

expect_s3_class(asr_cluster, "data.frame")

expect_equal(object = nrow(asr_cluster),expected = length(tr$tip.label))
