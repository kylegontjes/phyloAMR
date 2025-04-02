library(testthat)
library(phylosuite)
library(dplyr)
library(stringr)
library(ggplot2)
library(ape)
tr <- phylosuite::tr
df <- phylosuite::df
df$isolate_no <- df$tip_name_var
pheno <-  'colistin_ns'
tip_name_var <- "tip_name_var"

asr_obj <- asr(df = df,tr = tr,tip_name_var = tip_name_var ,pheno = pheno,model="ARD",node_states = "joint")
asr_cluster <- asr_cluster_detection(df = df,tr = tr,pheno = "colistin_ns",tip_name_var = "tip_name_var",patient_id = "PatientID",node_states = 'joint',faux_clusters = 'remove',parent_child_df = asr_obj$parent_child_df,confidence = NULL,remove_revertant = "yes",collapse_cluster = "yes")

expect_s3_class(asr_cluster, "data.frame")

expect_equal(object = nrow(asr_cluster),expected = length(tr$tip.label))
