library(testthat)
library(tidyverse)
tr <- phylosuite::tr
df <- phylosuite::df
pheno <-  'colistin_ns'
tip_name_var <- "tip_name_var"

best_model <- phylosuite::find_best_model(df,tr,tip_name_var,pheno,"joint")

expect_type(best_model, "character")
expect(best_model %in% c("ER","ARD"),"Not in desired model list")
