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
#' @return Synchronous gain and loss events of two traits with p-value permutation testing results
#' @export
synchronous_permutation_test <- function(comparitor, df, tr, tip_name_var, trait, node_states = "joint", conf_threshold = NULL, confidence = NULL, num_permutations = 1000, num_cores = 6) {
  library(pbmcapply)
  # Trait
  trait_asr <- asr(df = df, tr = tr, tip_name_var = tip_name_var, pheno = pheno, model = "ER", node_states = node_states, conf_threshold = conf_threshold) %>% .$parent_child_df

  # Comparitor
  comparitor_asr <- asr(df = df, tr = tr, tip_name_var = tip_name_var, pheno = comparitor, model = "ER", node_states = node_states, conf_threshold = conf_threshold) %>% .$parent_child_df

  # Observed data
  trait_synchronous_obs <- synchronous_detection(comparitor_asr, trait_asr, node_states = node_states, confidence = confidence)
  # Expected data from permutatino testing
  num_isolates <- nrow(df)
  num_features <- sum(as.numeric(df[[genotype]]))
  trait_runs <- pbmclapply(seq_len(num_permutations), FUN = function(x) {
    sample(df[[genotype]], num_isolates, replace = FALSE)
  }, mc.cores = num_cores) %>% do.call(cbind, .) %>% data.frame
  rownames(trait_runs) <- df[, tip_name_var]
  permutation_names <- paste0("p", seq_len(num_permutations))
  colnames(trait_runs) <- permutation_names
  trait_runs <- data.frame(trait_runs, tip_name_var = df[, tip_name_var])
  asr_permutation <- pbmclapply(permutation_names, FUN = function(x) {
    asr_result <- asr(df = trait_runs, tr = tr, tip_name_var = "tip_name_var", pheno = x, model = "ER", node_states = "joint") %>% .$parent_child_df
    comparitor_sychronous_perm <- synchronous_detection(asr_result, pheno_asr)
    return(comparitor_sychronous_perm)
  }, mc.cores = num_cores) %>% do.call(rbind, .)

  # Calculate p-values for each event of interest and add them to the trait observation data
  trait_synchronous_obs$synchronous_gain_pval <- c(1 + sum(asr_permutation$synchronous_gains_num >= trait_synchronous_obs$synchronous_gains_num)) / c(1 + num_permutations)
  trait_synchronous_obs$synchronous_loss_pval <- c(1 + sum(asr_permutation$synchronous_loss_num >= trait_synchronous_obs$synchronous_losses_num)) / c(1 + num_permutations)
  trait_synchronous_obs$synchronous_gain_loss_pval <- c(1 + sum(asr_permutation$synchronous_gain_loss_num >= trait_synchronous_obs$synchronous_gain_loss_num)) / c(1 + num_permutations)
  trait_synchronous_obs$synchronous_loss_gain_pval <- c(1 + sum(asr_permutation$synchronous_loss_gain_num >= trait_synchronous_obs$synchronous_loss_gain_num)) / c(1 + num_permutations)

  return(trait_synchronous_obs)
}
