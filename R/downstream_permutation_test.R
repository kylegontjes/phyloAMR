#' downstream_permutation_test: Permutation test for downstream analysis
#'
#' Permutation test for analysis of downstream gain and loss of a comparitor trait on stretches of a different trait of interest
#'
#' @param comparitor Comparitor trait, such as a genotype or a different trait/phenotype
#' @param df Dataframe with trait and comparitor and tip_name_variable
#' @param tr Phylogenetic tree
#' @param tip_name_variable Tip name variable
#' @param trait Trait of interest. Their stretches of trait presence will be characterized and evaluated
#' @param node_states Joint or marginal reconstruction
#' @param confidence_threshold Set value for confidence threshold of MLE calls if marginal
#' @param confidence Whether to use high or low confidence transition nodes when node_states are marginal
#' @param num_permutations Number of permutations. Default = 1000
#' @param num_cores Number of cores. Default = 6
#' @return Summary stats for downstream gain and loss of a trait
#' @export
downstream_permutation_test <- function(comparitor, df, tr, tip_name_variable, trait, node_states = "joint", confidence_threshold = NULL, confidence = NULL, num_permutations = 1000, num_cores = 6) {
  tr$node.label <- NULL

  # trait
  trait_asr <- asr(df = df, tr = tr, tip_name_variable = tip_name_variable, trait = trait, model = "ER", node_states = node_states, confidence_threshold = confidence_threshold)
  trait_asr_parent_child <- trait_asr$parent_child_df

  # Comparitor
  comparitor_asr <- asr(df = df, tr = tr, tip_name_variable = tip_name_variable, trait = comparitor, model = "ER", node_states = node_states, confidence_threshold = confidence_threshold)
  comparitor_asr_parent_child <- comparitor_asr$parent_child_df

  # Get observed
  downstream_nodes <- get_trait_traces_on_tree(parent_child_df = trait_asr_parent_child, tr = tr, node_states = node_states)
  downstream_changes <- get_gain_loss_on_stretches(comparitor_parent_child_df = comparitor_asr_parent_child, downstream_nodes = downstream_nodes, node_states = node_states, confidence = confidence)

  # Break if transitions = 0
  if (downstream_changes$transitions_num == 0) {
    downstream_changes$transition_pval <-  1
    downstream_changes$transition_stretches_pval <- 1
    downstream_changes$gain_pval <- 1
    downstream_changes$gains_stretches_pval <- 1
    downstream_changes$loss_pval <- 1
    downstream_changes$loss_stretches_pval <- 1
    asr_permutation <- NULL
  } else {

  # Generate rate data
  comparitor_index_mat <- comparitor_asr$corHMM_out$index.mat
  comparitor_p <- sapply(1:max(comparitor_index_mat, na.rm = TRUE), function(x)
      na.omit(c(comparitor_asr$corHMM_out$solution))[na.omit(c(comparitor_index_mat) == x)][1])

  # Expected data from permutation testing
  num_isolates <- nrow(df)
  comparitor_vals <- df[[comparitor]]
  trait_runs <- replicate(num_permutations, comparitor_vals[sample.int(num_isolates)], simplify = FALSE)

  # Permutation test
  asr_permutation <- parallel::mclapply(trait_runs, FUN = function(x) {
    dataset <- cbind(tip_names, x)
    outcome_str <- setNames(x, tip_names)
    asr_recon <- corHMM::ancRECON(tr, dataset, p = comparitor_p, method = node_states, rate.cat = 1, rate.mat = comparitor_asr$corHMM_out$index.mat, root.p = comparitor_asr$corHMM_out$root.p, get.likelihood = FALSE, get.tip.states = FALSE)
    asr_result <- get_parent_child_data(tr = tr, ancestral_states = asr_recon$lik.anc.states, trait_data = outcome_str, confidence_threshold = confidence_threshold, node_states = node_states)
    asr_parent_child_df <- get_continuation_data(asr_result, node_states)
    downstream_perm <- get_gain_loss_on_stretches(comparitor_parent_child_df = asr_parent_child_df, downstream_nodes = downstream_nodes, node_states = node_states, confidence = confidence)
    return(downstream_perm)
  }, mc.cores = num_cores)
  asr_permutation <- do.call(rbind, asr_permutation)

  downstream_changes$transition_pval <-  c(1 + sum(asr_permutation$transitions_num >= downstream_changes$transitions_num)) / c(1 + num_permutations)
  downstream_changes$transition_stretches_pval <-  c(1 + sum(asr_permutation$stretches_w_transitions_num >= downstream_changes$stretches_w_transitions_num)) / c(1 + num_permutations)
  downstream_changes$gain_pval <- c(1 + sum(asr_permutation$gains_num >= downstream_changes$gains_num)) / c(1 + num_permutations)
  downstream_changes$gains_stretches_pval <- c(1 + sum(asr_permutation$stretches_w_gains_num >= downstream_changes$stretches_w_gains_num)) / c(1 + num_permutations)
  downstream_changes$loss_pval <- c(1 + sum(asr_permutation$loss_num >= downstream_changes$loss_num)) / c(1 + num_permutations)
  downstream_changes$loss_stretches_pval <- c(1 + sum(asr_permutation$stretches_w_losses_num >= downstream_changes$stretches_w_losses_num)) / c(1 + num_permutations)
  }

  results <- list(downstream_permutation_testing = downstream_changes, permutations = asr_permutation)

  return(results)
}
