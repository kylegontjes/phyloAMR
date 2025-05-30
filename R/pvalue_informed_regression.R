#' Logistic regression using p-values to select final model
#'
#' Performs unadjusted logistic regression to identify candidate variables that fall under an p-value threshold (entry_criteria). Forward variable selection is performed to introduce variables into the model and retain if they fall within a more stringent criteria (retention_criteria).
#'
#' @param outcome Outcome of interest
#' @param variables Exposure variables of interest. Must be numeric or one-hot encoded
#' @param dataset Dataframe that contains the trait and exposure variables
#' @param entry_criteria P-value criteria for entry into the model. Default = 0.2
#' @param retention_criteria  P-value criteria for retention into the model. Default = 0.1
#' @return Final regression model
#' @importFrom stats terms
pvalue_informed_regression <- function(outcome, dataset, variables, entry_criteria = 0.2, retention_criteria = 0.1) {

  # Step 1: Univariable regression
  univariable <- univariable_regression(outcome = outcome, variables = variables, dataset = dataset)

  # Step 2: Identify variables under entry threshold
  univariable <- arrange(univariable, p_value)
  pvalue_under_threshold <- univariable[univariable$p_value < as.numeric(entry_criteria), ]

  # Step 3: Loop through variables to create final model
  ## Check if at least one variable is under criteria
  if(nrow(pvalue_under_threshold) ==0) {
    final_model <- paste0("No variables under entry criteria: ", entry_criteria)
  } else {
  # If at least one variable is under criteria generate model with
    eligible_variables <- pvalue_under_threshold$variable
    # Now loop through variables
    model <- paste(outcome, "1", sep = "~")
    glm_model <- glm(model, data = dataset, family = "binomial")
    for(i in eligible_variables){
      model_new <- paste0(glm_model$formula, "+", i)
      glm_model_new <- glm(model_new, data = dataset, family = "binomial")
      pvals <- summary(glm_model_new)$coefficients[, "Pr(>|z|)"]
      if(pvals[length(pvals)] < retention_criteria) {
        glm_model <- glm_model_new
      }
    }

    # Step 4: Format regression model
    ## Check if final model has any variables
    if(length(attr(terms(glm_model), "term.labels")) == 0){
      final_model <- paste0("No variables under retention criteria: ", retention_criteria)
    } else {
      # Format regression model
      final_model <- format_logistic_regression_table(glm_model)
    }

    # Step 5: Return results
    results <- list(final_model = final_model, univariable = univariable)
  }

  return(results)
}
