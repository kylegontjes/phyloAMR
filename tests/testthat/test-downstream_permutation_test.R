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
geno <-  colnames(genotype_mat)[3]
tip_name_variable <- "tip_name_var"

downstream_results <- downstream_permutation_test(comparitor = geno, df = df,tr =  tr, tip_name_variable = tip_name_variable,trait = trait, node_states = 'joint',confidence_threshold = NULL,num_permutations = 10,confidence = NULL, num_cores = 1)

expect_equal(object = nrow(downstream_results$downstream_permutation_testing),expected = 1)
