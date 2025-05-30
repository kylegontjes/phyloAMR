#' Logistic regression using stepwise Akaike information criteria selection
#'
#' Performs Akaike information criterion (AIC) informed variable selection using the stats::step function.
#'
#' @param outcome Outcome of interest. Character string.
#' @param dataset Dataframe that contains the trait and exposure variables
#' @param variables Exposure variables of interest. Must be numeric or one-hot encoded.
#' @param stepwise_direction  Direction of stepwise selection. Options include: 'both', 'backward', or 'forward'. For more information on these selection methods, see stats::step.
#' @return Regression results
#' @importFrom stats terms
AIC_stepwise_regression <- function(outcome, dataset, variables, stepwise_direction) {
  # Full model
  model <- as.formula(paste0(outcome, " ~ 1 +", paste0(variables, collapse = " + ")))

  # Selection
  if (stepwise_direction == "forward") {
    # Forward selection
    initialize_model <- glm(formula =  as.formula(paste0(outcome, " ~ 1")), data = dataset, family = "binomial")
    stepwise_glm_model <- suppressMessages(stats::step(object = initialize_model, direction = stepwise_direction, scope = model, trace = 0))
  }   else {
    # Backward or Both selection
    glm_model <- glm(model, data = dataset, family = "binomial")
    stepwise_glm_model <- suppressMessages(stats::step(object = glm_model, direction = stepwise_direction, trace = 0))
  }

  # Format and report multivariable table
  ## Check if final model has any variables
  if(length(attr(terms(stepwise_glm_model), "term.labels")) == 0){
    results <- paste0("No variables retained in final model")
  } else {
    # Format regression model
    results <- format_logistic_regression_table(stepwise_glm_model)
  }

  return(results)
}
