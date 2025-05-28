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
