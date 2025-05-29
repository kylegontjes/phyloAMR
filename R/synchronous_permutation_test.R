#' Synchronous permutation test
#'
#' Permutation test for synchronous detection function
#'
#' @param comparitor Comparitor, such as a genotype or a different trait/phenotype
#' @param df Dataframe with comparitor, trait, and tip_name_variable
#' @param tr Phylogenetic tree
#' @param tip_name_variable Tip name variable
#' @param trait Trait of interest, such as a trait/phenotype
#' @param node_states Joint or marginal reconstruction
#' @param confidence_threshold Confidence threshold if using marginal reconstruction
#' @param confidence Whether to use high or low confidence transition nodes when node_states are marginal
#' @param num_permutations Number of permutations. Default: 1000
#' @param num_cores Number of cores. Default: 1
#' @return Synchronous gain and loss events of two traits with p-value permutation testing results
#' @export
synchronous_permutation_test <- function(comparitor, df, tr, tip_name_variable, trait, node_states = "joint", confidence_threshold = NULL, confidence = NULL, num_permutations = 1000, num_cores = 1) {
  df <- df[match(tr$tip.label, df[[tip_name_variable]]), ]
  tr$node.label <- NULL

  # Trait
  trait_asr <- asr(df = df, tr = tr, tip_name_variable = tip_name_variable, trait = trait, model = "ER", node_states = node_states, confidence_threshold = confidence_threshold)
  trait_asr_parent_child <- trait_asr$parent_child_df

  # Comparitor
  comparitor_asr <- asr(df = df, tr = tr, tip_name_variable = tip_name_variable, trait = comparitor, model = "ER", node_states = node_states, confidence_threshold = confidence_threshold)
  comparitor_asr_parent_child <- comparitor_asr$parent_child_df

  # Observed data
  trait_synchronous_obs <- synchronous_transitions(comparitor_asr_parent_child, trait_asr_parent_child, node_states = node_states, confidence = confidence)

  if (trait_synchronous_obs$synchronous_transitions_num == 0) {
    trait_synchronous_obs$transition_pval <-  1
    trait_synchronous_obs$synchronous_gain_pval <-  1
    trait_synchronous_obs$synchronous_loss_pval <-  1
    trait_synchronous_obs$synchronous_gain_loss_pval <-  1
    trait_synchronous_obs$synchronous_loss_gain_pval <-  1
    asr_permutation <- NULL
  } else {

  # Generate rate data
  comparitor_index_mat <- comparitor_asr$corHMM_out$index.mat
  comparitor_p <- sapply(1:max(comparitor_index_mat, na.rm = TRUE), function(x) {
    na.omit(c(comparitor_asr$corHMM_out$solution))[na.omit(c(comparitor_index_mat) == x)][1]
    })

  # Expected data from permutation testing
  num_isolates <- nrow(df)
  comparitor_vals <- df[[comparitor]]
  trait_runs <- replicate(num_permutations, comparitor_vals[sample.int(num_isolates)], simplify = FALSE)
  tip_names <- df[[tip_name_variable]]

  # Permutations
  asr_permutation <- parallel::mclapply(trait_runs, FUN = function(x) {
    dataset <- cbind(tip_names, x)
    outcome_str <- setNames(x, tip_names)
    asr_recon <- corHMM::ancRECON(tr, dataset, p = comparitor_p, method = node_states, rate.cat = 1, rate.mat = comparitor_asr$corHMM_out$index.mat, root.p = comparitor_asr$corHMM_out$root.p, get.likelihood = FALSE, get.tip.states = FALSE)
    asr_parent_child <- get_parent_child_data(tr = tr, ancestral_states = asr_recon$lik.anc.states, trait_data = outcome_str, confidence_threshold = confidence_threshold, node_states = node_states)
    asr_result <- get_continuation_data(parent_child_df = asr_parent_child, node_states = node_states)
    comparitor_sychronous_perm <- synchronous_transitions(comparitor_parent_child_df = asr_result, trait_parent_child_df =  trait_asr_parent_child, node_states = node_states, confidence = confidence)
    return(comparitor_sychronous_perm)
  }, mc.cores = num_cores)
  asr_permutation <- do.call(rbind, asr_permutation)

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
