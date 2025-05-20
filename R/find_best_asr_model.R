#' find_best_asr_model: Determine best model for ancestral state reconstruction
#'
#' Function to find best model for ancestral state reconstruction using corHMM.
#'
#' @param df Dataframe with tip name (e.g., tip_name_variable) and phenotype/trait (e.g., trait) variables
#' @param tr Phylogenetic tree
#' @param tip_name_variable Name of variable containing tip names in df
#' @param trait Name of phenotype/trait variable in df
#' @param node_states Whether to perform "joint" or "marginal" reconstruction
#' @param upper_bound Upper bound for likelihood search. Default: 1e50
#' @param lower_bound Lower bound for likelihood search. Default: 1e-9
#' @return Description of return value
#'   \describe{
#'     \item{model_statistics}{Summary of model options}
#'     \item{best_model}{Best model for ancestral state reconstruction}
#'     \item{node_states}{Whether "joint" or "marginal" reconstruction was performed}
#'     \item{corHMM_output_best_model}{corHMM output for best model}
#'   }
#' @export

find_best_asr_model <- function(df, tr, tip_name_variable, trait, node_states = "joint", upper_bound = 1e50, lower_bound = 1e-9) {
  # Checks
  ## Check if trait is 0,1
  check_trait(df[[trait]])

  # Run corHMM to estimate hidden rates
  ## Equal rates (ER) model
  cat("Testing equal rates (ER) model \n")
  ### Generate corHMM model
  corHMM_ER <- invisible(corHMM::corHMM(phy = tr, data = df[, c(tip_name_variable, trait)], rate.cat = 1, model = "ER", node.states = node_states, upper.bound = upper_bound, lower.bound = lower_bound))
  ### Generate model statistics
  corHMM_ER_stats <- characterize_asr_model(corHMM_ER)
  cat(paste0("AICc of ER model is: ", corHMM_ER_stats$AICc, " \n"))

  ## All rates differ (ARD) model
  cat("Testing all rates differ (ARD) model \n")
  ### Generate corHMM model
  corHMM_ARD <- invisible(corHMM::corHMM(phy = tr, data = df[, c(tip_name_variable, trait)], rate.cat = 1, model = "ARD", node.states = node_states, upper.bound = upper_bound, lower.bound = lower_bound))
  ### Generate model statistics
  corHMM_ARD_stats <- characterize_asr_model(corHMM_ARD)
  cat(paste0("AICc of ARD model is: ", corHMM_ARD_stats$AICc, " \n"))

  # Model options table
  model_statistics <- rbind(corHMM_ER_stats, corHMM_ARD_stats)

  # Identify best model with the lowest AICc
  best_model  <- model_statistics[which.min(model_statistics$AICc), "model"]
  model_statistics$best_model <- ifelse(model_statistics$model == best_model, TRUE, FALSE)
  cat(paste0("Best model: ", best_model, " \n"))

  # Store best model's output
  corHMM_output_best_model <- if (best_model == "ER") {
    corHMM_ER
  } else {
    corHMM_ARD
  }

  # Results object
  results <- list(model_statistics = model_statistics,
                  best_model = best_model,
                  node_states = node_states,
                  corHMM_output_best_model = corHMM_output_best_model)

  return(results)
}
