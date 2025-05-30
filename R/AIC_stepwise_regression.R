#' Logistic regression using stepwise Akaike information criteria selection
#'
#' Performs Akaike information criterion (AIC) informed selection using stats::step function.
#'
#' @param outcome Outcome of interest
#' @param dataset Dataframe that contains the trait and exposure variables
#' @param variables Exposure variables of interest. Must be numeric or one-hot encoded
#' @param stepwise_direction  Direction of stepwise selection. Options include: 'both', 'backward', or 'forward'. For more information, see stats::step.
#' @return Regression results

AIC_stepwise_regression <- function(outcome, dataset, variables, stepwise_direction) {
  if (stepwise_direction == "forward") {
    # Forward selectoin
    initialize_model <- glm(formula =  as.formula(paste0(outcome, " ~ 1")), data = dataset, family = "binomial")
    stepwise_glm_model <- suppressMessages(stats::step(initialize_model, direction = stepwise_direction, scope = model, trace = 0))
  }   else {
    # Backward or Both selection
    model <- as.formula(paste0(outcome, " ~ 1 +", paste0(variables, collapse = " + ")))
    glm_model <- glm(model, data = dataset, family = "binomial")
    stepwise_glm_model <- suppressMessages(stats::step(glm_model, direction = stepwise_direction, trace = 0))
  }

  # Format results
  results <- format_logistic_regression_table(stepwise_glm_model)
  return(results)
}
