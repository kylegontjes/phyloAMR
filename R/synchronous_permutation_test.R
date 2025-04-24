#' synchronous_permutation_test: Synchronous permutation test
#'
#' Permutation test for synchronous detection function
#'
#' @param comparitor Comparitor, such as a genotype
#' @param df Dataframe with comparitor, trait, and tip_name_var
#' @param tr Phylogenetic tree
#' @param tip_name_var Tip name variable
#' @param trait Trait of interest, such as a phenotype
#' @param node_states Joint or marginal reconstruction
#' @param conf_threshold Confidence threshold if using marginal reconstruction
#' @param confidence Whether to use high or low confidence transition nodes when node_states are marginal
#' @param num_permutations Number of permutations. Default: 1000
#' @param num_cores Number of cores. Default: 1
#' @return Synchronous gain and loss events of two traits with p-value permutation testing results
#' @export
synchronous_permutation_test <- function(comparitor, df, tr, tip_name_var, trait, node_states = "joint", conf_threshold = NULL, confidence = NULL, num_permutations = 1000, num_cores = 1) {
  df <- df %>% .[match(tr$tip.label, .[[tip_name_var]]), ]
  tr$node.label <- NULL

  # Trait
  trait_asr <- asr(df = df, tr = tr, tip_name_var = tip_name_var, pheno = pheno, model = "ER", node_states = node_states, conf_threshold = conf_threshold) %>% .$parent_child_df

  # Comparitor
  comparitor_asr <- asr(df = df, tr = tr, tip_name_var = tip_name_var, pheno = comparitor, model = "ER", node_states = node_states, conf_threshold = conf_threshold)
  comparitor_asr_parent_child <- comparitor_asr$parent_child_df

  # Observed data
  trait_synchronous_obs <- synchronous_detection(comparitor_asr_parent_child, trait_asr, node_states = node_states, confidence = confidence)

  if(trait_synchronous_obs$synchronous_transitions_num == 0) {
    trait_synchronous_obs$transition_pval <-  1
    trait_synchronous_obs$synchronous_gain_pval <-  1
    trait_synchronous_obs$synchronous_loss_pval <-  1
    trait_synchronous_obs$synchronous_gain_loss_pval <-  1
    trait_synchronous_obs$synchronous_loss_gain_pval <-  1
    asr_permutation <- NULL
  } else {

    # Generate rate data
  comparitor_asr_corHMM <- comparitor_asr$corHMM_out
  comparitor_index_mat <- comparitor_asr$corHMM_out$index.mat
  comparitor_solution_mat <- comparitor_asr$corHMM_out$solution
  comparitor_root.p <- comparitor_asr$corHMM_out$root.p
  comparitor_p <- sapply(1:max(comparitor_index_mat, na.rm = TRUE), function(x)
      na.omit(c(comparitor_solution_mat))[na.omit(c(comparitor_index_mat) == x)][1])

  # Expected data from permutation testing
  num_isolates <- nrow(df)
  trait_runs <- pbmclapply(seq_len(num_permutations), FUN = function(x) {
    sample(df[[comparitor]], num_isolates, replace = FALSE)
  }, mc.cores = num_cores) %>% do.call(cbind, .) %>% data.frame
  permutation_names <- paste0("T", seq_len(num_permutations))
  colnames(trait_runs) <- permutation_names
  trait_runs <- data.frame(trait_runs %>% mutate_all(as.integer), tip_name_var = df[, tip_name_var])

  # Permutations
  asr_permutation <- pbmclapply(permutation_names, FUN = function(x) {
    dataset =  trait_runs[,c('tip_name_var',x)]
    outcome_str <- trait_runs[, x] %>% `names<-`(trait_runs[['tip_name_var']])
    asr_recon <- ancRECON(tr,dataset,p = comparitor_p,method = 'joint',rate.cat = 1,rate.mat = comparitor_asr_corHMM$index.mat)
    asr_result <- get_parent_child_data(tr = tr,anc_data = asr_recon$lik.anc.states,pheno_data = outcome_str,conf_threshold = NULL,node_states = 'joint') %>% get_continuation_data(.,'joint')
    comparitor_sychronous_perm <- synchronous_detection(asr_result, trait_asr)
    return(comparitor_sychronous_perm)
  }, mc.cores = num_cores) %>% do.call(rbind, .)

  # Calculate p-values for each event of interest and add them to the trait observation data
  trait_synchronous_obs$transition_pval <- c(1 + sum(asr_permutation$synchronous_transitions_num >= trait_synchronous_obs$synchronous_transitions_num)) / c(1 + num_permutations)
  trait_synchronous_obs$synchronous_gain_pval <- c(1 + sum(asr_permutation$synchronous_gains_num >= trait_synchronous_obs$synchronous_gains_num)) / c(1 + num_permutations)
  trait_synchronous_obs$synchronous_loss_pval <- c(1 + sum(asr_permutation$synchronous_losses_num >= trait_synchronous_obs$synchronous_losses_num)) / c(1 + num_permutations)
  trait_synchronous_obs$synchronous_gain_loss_pval <- c(1 + sum(asr_permutation$synchronous_gain_loss_num >= trait_synchronous_obs$synchronous_gain_loss_num)) / c(1 + num_permutations)
  trait_synchronous_obs$synchronous_loss_gain_pval <- c(1 + sum(asr_permutation$synchronous_loss_gain_num >= trait_synchronous_obs$synchronous_loss_gain_num)) / c(1 + num_permutations)

  }
  results <- list(synchronous_permutation_testing = trait_synchronous_obs, permutations = asr_permutation)
  return(results)
}
