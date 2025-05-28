pvalue_informed_regression <- function(outcome, dataset, variables, entry_criteria = 0.2, retention_criteria = 0.1) {
  # Univariable regression
  univariable <- univariable_regression(outcome = outcome, variables = variables, dataset = dataset)
  # Variables under entry threshold
  univariable <- arrange(univariable, p_value)
  pval_less_threshold <- datatable[univariable$p_value < as.numeric(entry_criteria), ]
  eligible_variables <- pval_less_threshold$variable

  # Loop through variables to create final model
  ## Check if at least one variable is under criteria
  if(min(pval_less_threshold$p_value) > retention_criteria) {
    return(paste0("No model - no variables < ", retention_criteria))
  } else {
  # If at least one variable is under criteria generate model with
    # Now loop through variables
    model <- paste(outcome, "1", sep = "~")
    glm_model <- glm(model, data = dataset, family = "binomial")
    # Retain if pvalue is
    best_variable <- pval_less_threshold[which.min(pval_less_threshold$p_value),"variable"]
    model <- paste0(glm_model$formula, "+", best_variable)
    glm_model <- glm(model, data = dataset, family = "binomial")
    for(i in eligible_variables[2:length(eligible_variables)]){
      model_new <- paste0(glm_model$formula, "+", i)
      glm_model_new <- glm(model_new, data = dataset, family = "binomial")
      pvals <- round(abs(summary(glm_model_new)$coefficients[, "Pr(>|z|)"]), 4)
      if(pvals[length(pvals)] > retention_criteria) {
        glm_model <- glm_model
      } else {
        glm_model <- glm_model_new
      }
    }

    # Format regression model
    final_model <- format_logistic_regression_table(glm_model)

    # Results
    results <- list(final_model = final_model, univariable = univariable)
    return(results)
  }
}
