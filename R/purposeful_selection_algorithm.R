#' Logistic regression using purposeful selection
#'
#' Uses purposeful selection algorithm to identify a regression model.
#' Three steps exist:
#' 1. Unadjusted logistic regression to identify candidate variables under a p-value threshold (entry_criteria)
#' 2. Multivariable regression of candidate variables. Iterative process, starting from the value with the highest p-value, variables are retained if they fall under category of a significant variable (< retention_criteria) or confounding (i.e., effect size +/- confounding criteria).
#' 3. All variables that failed step 1 are introduced to the model and retained if p-value < retention criteria (retention_criteria)
#' Reference: https://scfbm.biomedcentral.com/articles/10.1186/1751-0473-3-17
#'
#' @param outcome Outcome of interest
#' @param variables Exposure variables of interest. Must be numeric or one-hot encoded
#' @param dataset Dataframe that contains the trait and exposure variables
#' @param entry_criteria P-value criteria for entry into the model. Default = 0.2
#' @param retention_criteria  P-value criteria for retention into the model. Default = 0.1
#' @param confounding_criteria  Percent change of effect size. Set value to be very high (i.e., 1000) if testing for confounding is not desired. Default = 0.2.
#' @return A list containing univariable results (ps_step_1), initial multivariable modeling (ps_step_2), model refinement with non-candidate variables (ps_step_3), and the final glm model (final_model) and table (final_model_table).
#' @importFrom dplyr left_join
#' @importFrom stats formula
#' @importFrom stats terms
purposeful_selection_algorithm <- function(outcome, variables, dataset, entry_criteria = 0.2, retention_criteria = 0.1, confounding_criteria = 0.2) {
  # Step 1. Unadjusted regression to identify candidate variables
  ps_step_1 <- purposeful_selection_step_1(outcome = outcome, variables = variables, dataset = dataset, entry_criteria = entry_criteria)

  # Step 2. Multivariable regression of candidate variables and iterative variable retention analysis
  ps_step_2 <- purposeful_selection_step_2(outcome = outcome, candidate_variables = ps_step_1$candidates, dataset = dataset, retention_criteria = retention_criteria, confounding_criteria = confounding_criteria)

  # Step 3. Evaluation of non-candidate variables from step 1 in final model from step 2.
  if (length(ps_step_2) == 1 && sum(ps_step_2 == "No significant variables") > 0) {
    results <- list(ps_step_1 = ps_step_1, warning = "No multivariable model possible as all candidate variables were not statistically significant.")
  } else {
    p3_candidate_variables <- variables[!variables %in% ps_step_1$candidates]
    ps_step_3 <- purposeful_selection_step_3(outcome = outcome, fixed_model_variables = ps_step_2$model_variables, candidate_variables = p3_candidate_variables, dataset = dataset, retention_criteria = retention_criteria, confounding_criteria = confounding_criteria)
    results <- list(ps_step_1 = ps_step_1, ps_step_2 = ps_step_2, ps_step_3 = ps_step_3, final_model = ps_step_3$final_model, final_model_table = ps_step_3$final_model_table)
  }
  return(results)
}

purposeful_selection_step_1 <- function(outcome, variables, dataset, entry_criteria) {
  # Univariable regression
  univariable_results <- univariable_regression(outcome = outcome, dataset = dataset, variables = variables)

  # Identify candidate and non-candidate variables
  candidates <- univariable_results[univariable_results[["p_value"]] < entry_criteria, "variable"]
  non_candidates <- subset(variables, !variables %in% candidates)

  # Return regression results and variable lists
  results <- list(candidates = candidates, non_candidates = non_candidates, univariable_regression = univariable_results)
  return(results)
}

purposeful_selection_step_2 <- function(outcome, candidate_variables, dataset, retention_criteria, confounding_criteria) {
  # Initial multivariable model with all variables included
  input_formula <- paste0(outcome, " ~ 1 + ", paste0(candidate_variables, collapse = " + "))
  glm_model <- glm(formula = input_formula, data = dataset, family = "binomial")

  # Iterative backward selection on candidate variables to get initial
  iterative_modeling <- iterative_significance_model_testing(glm_model = glm_model, candidate_variables = candidate_variables, fixed_variables = NULL, dataset = dataset, retention_criteria = retention_criteria, confounding_criteria =  confounding_criteria)

  return(iterative_modeling)
}

purposeful_selection_step_3 <-  function(outcome, fixed_model_variables, candidate_variables, dataset, retention_criteria, confounding_criteria) {
  # Input formula
  input_formula <- paste0(outcome, " ~ 1 + ", paste0(fixed_model_variables, collapse = " + "))
  # Determine candidates for retention analysis
  model_additions <- lapply(candidate_variables, FUN = function(x) {
    glm_model <- glm(formula = paste0(input_formula, " + ", x), data = dataset, family = "binomial")
    format_logistic_regression_table(glm_model)
  })
  # Curate regression results for untested variables
  model_additions_table <- do.call(rbind, model_additions)
  model_additions_table <- subset(model_additions_table, !variable %in% c("(Intercept)", fixed_model_variables))
  model_additions_table$p_value <- as.numeric(model_additions_table$p_value)
  model_additions_table <- arrange(model_additions_table, -p_value)
  # Final candidates
  final_candidates <- model_additions_table[model_additions_table[["p_value"]] < retention_criteria, "variable"]
  if (length(final_candidates) == 0) {
    glm_model <- glm(formula = input_formula, data = dataset, family = "binomial")
    glm_model_table <- format_logistic_regression_table(glm_model)
    # Results
    results <- list(candidate_variables = candidate_variables, model_additions_table = model_additions_table, final_model = glm_model, final_model_table = glm_model_table)
  } else {
    updated_formula <- paste0(input_formula, "+", paste0(final_candidates, collapse = "+"))
    glm_model <- glm(formula = updated_formula, data = dataset, family = "binomial")

    # Iterative backward selection on candidate variables to get initial
    iterative_modeling <- iterative_significance_model_testing(glm_model = glm_model, candidate_variables = final_candidates, fixed_variables = fixed_model_variables, dataset = dataset, retention_criteria = retention_criteria, confounding_criteria =  confounding_criteria)

    # Results
    results <- list(candidate_variables = candidate_variables, model_additions_table = model_additions_table, final_model = iterative_modeling$glm_model, final_model_table = iterative_modeling$glm_model_table)
  }

  return(results)
}

iterative_significance_model_testing <- function(glm_model, candidate_variables, fixed_variables, dataset, retention_criteria, confounding_criteria) {
  # Generate glm model
  glm_model_table <- format_logistic_regression_table(glm_model)

  # Variables tested
  if (is.null(fixed_variables)) {
    tested <- c()
  } else {
    tested <- fixed_variables
  }

  confounder_internal <- c()

  # Iterative process: For top model component, drop if alpha lvel is above 0.1 AND not a confounder (two checks)
  repeat {
    # Determine if all candidate variables were tested
    if (sum(tested %in% candidate_variables) == length(candidate_variables)) {
      message("Tested all candidate variables")
      break
    }

    # Grab the p-value and name of variable with highest p-value that was not already tested
    glm_model_table_without_tested <- subset(glm_model_table, !variable %in% c(confounder_internal, tested))
    max_pval <- glm_model_table_without_tested[which.max(glm_model_table_without_tested[["p_value"]]), ]

    # If variable with max p-value is less than retention criteria, the iterative process concludes as model is compiled of either all
    if (max_pval[["p_value"]] < retention_criteria) {
      message("Retained variables met retention criteria or are confounders")
      # Model is completed as all either confounders or significant
      break
    } else {
      if (c(length(tested) + 1) == length(candidate_variables) && max_pval[["p_value"]] > retention_criteria) {
        warning("No significant variables, reporting step 1 analysis")
        results <- "No significant variables"
        break
      } else {
        # Test for confounding
        confounding <- test_confounding(variable = max_pval$variable, glm_model = glm_model, dataset = dataset, confounding_criteria = confounding_criteria)
        if (confounding$confounder_status == FALSE) {
          # Remove variable from consideration
          tested <- c(tested, max_pval$variable)
          confounder_internal <- c()
          # Remove variable from model
          glm_model <- confounding$reduced_model
          # Regenerate summary
          glm_model_table <- format_logistic_regression_table(glm_model)
        } else {
          glm_model <- glm_model
          confounder_internal <- c(confounder_internal, max_pval$variable)
        }
      }
    }
  }
  model_variables <- glm_model_table$variable
  results <- list(glm_model = glm_model, glm_model_table =  glm_model_table, model_variables = model_variables, fixed_variables = fixed_variables, confounder = confounder_internal)
  return(results)
}

test_confounding <- function(variable, glm_model, dataset, confounding_criteria) {
  # Remove term from model
  outcome <- all.vars(formula(glm_model))[[1]]
  model_terms <- attr(terms(glm_model), "term.labels")
  model_terms_without_variable <- model_terms[model_terms != variable]
  reduced_model_formula <-  paste0(outcome, " ~ 1 + ", paste0(model_terms_without_variable, collapse = " + "))
  # Reduced model
  glm_reduced_model <- glm(formula = reduced_model_formula, data = dataset, family = "binomial")
  # Generate variable with
  full_model_summary <-  format_logistic_regression_table(glm_model)[, c("variable", "OR")]
  colnames(full_model_summary) <- c("variable", "OR_full")
  reduced_model_summary <- format_logistic_regression_table(glm_reduced_model)[, c("variable", "OR")]
  colnames(reduced_model_summary) <- c("variable", "OR_reduced")
  # Compare models
  model_comparison <- suppressMessages(left_join(full_model_summary, reduced_model_summary))
  model_comparison <- model_comparison[model_comparison[["variable"]] != variable, ]
  # Convert OR to numeric
  model_comparison[["OR_full"]] <- as.numeric(model_comparison[["OR_full"]])
  model_comparison[["OR_reduced"]] <- as.numeric(model_comparison[["OR_reduced"]])
  # Test for confounding
  model_comparison$effect_change <- (model_comparison[["OR_reduced"]] - model_comparison[["OR_full"]]) / model_comparison[["OR_full"]]
  model_comparison$confounder <- ifelse(abs(model_comparison[["effect_change"]]) > confounding_criteria, "yes", "no")
  # IS it a confounder
  confounder_status <- ifelse(sum(model_comparison$confounder == "yes") > 0, TRUE, FALSE)
  # Results
  results <- list(confounder_status = confounder_status, confounder_table = model_comparison, reduced_model = glm_reduced_model)
  return(results)
}
