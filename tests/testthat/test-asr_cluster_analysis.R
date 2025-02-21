library(testthat)
library(tidyverse)
library(ape)
library(phytools)
tr <- phylosuite::tr
df <- phylosuite::df
df$isolate_no <- df$tip_name_var
pheno <-  'colistin_ns'
tip_name_var <- "tip_name_var"

asr_obj <- asr(df,tr,tip_name_var ,pheno,model="ARD",node_states = "joint",conf_threshold=0.875)
asr_cluster <- asr_cluster_detection(df = df,tr = tr,pheno = "colistin_ns",parent_child_df = asr_obj$parent_child_df,remove_faux = "yes",confidence = "high",remove_revertant = "yes",collapse_cluster = "yes")
asr_cluster_analysis <- asr_cluster_analysis(asr_cluster,"yes")

expect_s3_class(asr_cluster_analysis, "data.frame")

expect_equal(object = nrow(asr_cluster_analysis),expected = 1)
expect_equal(object = ncol(asr_cluster_analysis),expected = 27)
