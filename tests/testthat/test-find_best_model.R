library(testthat)
library(phylosuite)
library(dplyr)
library(stringr)
library(ggplot2)
tr <- phylosuite::tr
df <- phylosuite::df
pheno <-  'colistin_ns'
tip_name_var <- "tip_name_var"

best_model <- phylosuite::find_best_asr_model(df = df,tr = tr,tip_name_var = tip_name_var,pheno = pheno,node_states = "joint")

expect_type(best_model$best_model, "character")
expect(best_model$best_model %in% c("ER","ARD","SYS"),"Not in desired model list")
