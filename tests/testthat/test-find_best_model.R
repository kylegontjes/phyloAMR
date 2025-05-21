library(testthat)
library(phyloAMR)
library(dplyr)
library(stringr)
library(ggplot2)
tr <- phyloAMR::tr
df <- phyloAMR::df
trait <-  'colistin_ns'
tip_name_variable <- "tip_name_var"

best_model <- phyloAMR::find_best_asr_model(df = df, tr = tr, tip_name_variable = tip_name_variable, trait = trait, node_states = "joint")

expect_type(best_model$best_model, "character")
expect(best_model$best_model %in% c("ER", "ARD", "SYS"), "Not in desired model list")
