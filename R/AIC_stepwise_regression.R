#' Logistic regression using stepwise Akaike information criteria selection
#'
#' Performs Akaike information criteria (AIC) informed selection using stats::step function.
#'
#' @param outcome Outcome of interest
#' @param dataset Dataset with outcome and variables
#' @param variables Variables to evaluate. Must binarize variables (i.e., 0, 1).
#' @param stepwise_direction  Stepwise directions include "both", "backward", or "forward". For more information, see stats::step.
#' @return Regression results

AIC_stepwise_regression <- function(outcome, dataset, variables, stepwise_direction) {
  model <- as.formula(paste0(outcome, " ~ 1 +", paste0(variables, collapse = " + ")))
  if (stepwise_direction == "forward") {
    initialize_model <- glm(formula =  as.formula(paste0(outcome, " ~ 1")), data = dataset, family = "binomial")
    stepwise_glm_model <- suppressMessages(stats::step(initialize_model, direction = stepwise_direction, scope = model, trace = 0))
  }   else {
    glm_model <- glm(model, data = dataset, family = "binomial")
    stepwise_glm_model <- suppressMessages(stats::step(glm_model, direction = stepwise_direction, trace = 0))
  }
  results <- format_logistic_regression_table(stepwise_glm_model)
  return(results)
}
