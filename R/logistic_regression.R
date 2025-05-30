# Logistic regression base function regression
logistic_regression <- function(variables, outcome, dataset) {
  # Multvariable model
  model <- as.formula(paste0(outcome, " ~ 1 +", paste0(variables, collapse = " + ")))
  glm_model <- glm(formula = model, data = dataset, family = "binomial")

  # Format multivariable model
  results <- format_logistic_regression_table(glm_model)
  return(results)
}

# Univariable regression across numerous variables
univariable_regression <- function (outcome, dataset, variables) {
  datatable <- lapply(variables,
    FUN = function(x) {
      datatable <- logistic_regression(variables = x, outcome =  outcome, dataset =  dataset)
      return(datatable)
    }
  )
  results <- do.call(rbind.data.frame, datatable)
  return(results)
}

# Multivariable regression
multivariable_regression <- function(outcome, dataset, variables) {
  results <- logistic_regression(variables = variables, outcome = outcome, dataset = dataset)
  return(results)
}

format_logistic_regression_table <- function(glm_model) {
  # Coefficients
  ci <- suppressMessages(confint(glm_model))
  coefficients <- round(exp(cbind(OR = coef(glm_model), ci)), 2)
  coefficients <- format(coefficients, format = "f", digits = 4)

  # Pvalue
  p_value <- round(abs(summary(glm_model)$coefficients[, "Pr(>|z|)"]), 4)

  # Get variable name
  variable <-  all.vars(formula(glm_model))[-1]

  # Curate table with effect size and p-value. Also, remove the intercept results, cause not necessary!
  effect_size_pvalue <- cbind.data.frame(coefficients, p_value)
  effect_size_pvalue <- subset(effect_size_pvalue, rownames(effect_size_pvalue) != "(Intercept)")
  table <- cbind.data.frame(variable, effect_size_pvalue)
  colnames(table) <- c("variable","OR", "2.5%", "97.5%", "p_value")
  rownames(table) <- NULL

  # Arrange by p-value
  table <- arrange(table, p_value)
  table$p_value <- format(table$p_value, format = "f", digits = 4)

  return(table)
}
