#' find_best_asr_model: Determine best model for ancestral state reconstruction
#'
#' Function to find best model for ancestral state reconstruction using corHMM.
#'
#' @param df Dataframe with tip name variable and phenotype
#' @param tr Phylogenetic tree
#' @param tip_name_var Name of variable containing tip names in df
#' @param pheno Name of phenotype variable in df
#' @param node_states Whether to perform "joint" or "marginal" reconstruction
#' @param upper_bound Upper bound for likelihood search. Default: 1e50
#' @param lower_bound Lower bound for likelihood search. Default: 1e-9
#' @return Description of return value
#'   \describe{
#'     \item{model_options}{Summary of model options}
#'     \item{best_model}{Best model for ancestral state reconstruction}
#'   }
#' @export

find_best_asr_model <- function(df, tr, tip_name_var, pheno, node_states = "joint", upper_bound = 1e50, lower_bound = 1e-9) {
  # Check if phenotype is 0,1
  check_phenotype(df[[pheno]])

  # Run corHMM to estimate hidden rates
  corHMM_ER <- invisible(corHMM::corHMM(phy = tr, data = df[, c(tip_name_var, pheno)], rate.cat = 1, model = "ER", node.states = node_states, upper.bound = upper_bound, lower.bound = lower_bound) %>% characterize_asr_model())
  corHMM_ARD <- invisible(corHMM::corHMM(phy = tr, data = df[, c(tip_name_var, pheno)], rate.cat = 1, model = "ARD", node.states = node_states, upper.bound = upper_bound, lower.bound = lower_bound) %>% characterize_asr_model())

  model_options_tbl <- rbind(corHMM_ER, corHMM_ARD)

  best_model  <- model_options_tbl[which.min(model_options_tbl$AIC), "model"]

  results <- list(model_options = model_options_tbl,
                  best_model = best_model)
  return(results)
}
