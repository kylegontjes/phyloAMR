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
#' @param variables Variables to evaluate. Must binarize variables (i.e., 0, 1).
#' @param dataset Dataset with outcome and variables
#' @param entry_criteria P-value criteria for entry into the model. Default = 0.2
#' @param retention_criteria  P-value criteria for retention into the model. Default = 0.1
#' @param confounding_criteria  Percent change of effect size. Set value to be very high (i.e., 1000) if testing for confounding is not desired. Default = 0.2.
#' @return Regression results from each step.
#' @importFrom dplyr left_join
purposeful_selection_algorithm <- function(outcome, variables, dataset, entry_criteria = 0.2, retention_criteria = 0.1, confounding_criteria = 0.2) {
  ps_step1 <- purposeful_selection_step_1(outcome = outcome, variables = variables, dataset = dataset, entry_criteria = entry_criteria)
  ps_step2 <- purposeful_selection_step_2(outcome = outcome, candidate_variables = ps_step1$candidates, dataset = dataset, retention_criteria = retention_criteria, confounding_criteria = confounding_criteria)
  if(length(ps_step2) == 1 && sum(ps_step2 == "No significant variables") > 0) {
    results <- list(ps_step1 = ps_step1, warning = "No multivariable model possible as all candidate variables were not statistically significant.")
  } else {
    p3_candidate_variables <- subset(variables, !variables %in% ps_step1$candidates)
    ps_step3 <- purposeful_selection_step_3(outcome = outcome, fixed_model_variables = ps_step2$model_variables, candidate_variables = p3_candidate_variables, dataset = dataset, retention_criteria = retention_criteria, confounding_criteria = confounding_criteria)
    final_model_table <- format_logistic_regression_table(ps_step3$final_model)
    results <- list(ps_step1 = ps_step1, ps_step2 = ps_step2, ps_step3 = ps_step3, final_model = ps_step3$final_model, final_model_table = final_model_table)
  }
  return(results)
}

purposeful_selection_step_1 <- function(outcome, variables, dataset, entry_criteria) {
  univariable_results <- univariable_regression(outcome = outcome, dataset = dataset, variables = variables)
  candidates <- univariable_results[univariable_results[["p_value"]] < entry_criteria, "variable"]
  non_candidates <- subset(variables, !variables %in% candidates)
  results <- list(candidates = candidates, non_candidates = non_candidates, univariable_regression = univariable_results)
  return(results)
}

purposeful_selection_step_2 <- function(outcome, candidate_variables, dataset, retention_criteria, confounding_criteria) {
  input_formula <- paste0(outcome, " ~ 1 + ", paste0(candidate_variables, collapse = "+"))
  input_formula <- trimws(input_formula, whitespace = "\\+")
  # Glm model output to manipulate
  glm_model <- glm(formula = input_formula, data = dataset, family = "binomial")
  glm_model_tbl <- format_logistic_regression_table(glm_model)
  # Variables tested
  tested <- c()
  confounder <- c()
  # Iterative process: For top model component, drop if alpha lvel is above 0.1 AND not a confounder (two checks)
  repeat {
    # grab significant value
    max_pval <- glm_model_tbl[which.max(glm_model_tbl[["p_value"]]), ]
    if (max_pval[["p_value"]] < retention_criteria) {
      message("All candidate variables in step 2 met retention criteria")
      # Model is completed as all either confounders or significant
      break
    } else {
      if (c(length(tested) + 1) == length(candidate_variables) && max_pval[["p_value"]] > retention_criteria) {
        warning("No significant variables, reporting step 1 analysis")
        results <- "No significant variables"
        break
      } else {
        # Test as confounder
        ## Remove variable from model
        reduced_model_formula <- gsub(max_pval$variable, "", glm_model$formula)
        reduced_model_formula <- gsub("\\+\\+", "\\+", reduced_model_formula)
        reduced_model_formula <- trimws(reduced_model_formula, whitespace = "\\+")
        ## Reduced model
        glm_reduced_model <- glm(formula = reduced_model_formula, data = dataset, family = "binomial")
        # Test for confounding
        confounding <- test_confounding(variable = max_pval$variable, glm_model = glm_model, dataset = dataset, confounding_criteria = confounding_criteria)
        is_confounding <- ifelse(sum(confounding$confounder == "yes") > 0, "yes", "no")
        if (is_confounding == "no") {
          # Remove variable from model
          glm_model <- glm_reduced_model
          # Regenerate summary
          glm_model_tbl <- format_logistic_regression_table(glm_model)
          # Label tested variables
          tested <- c(tested, max_pval$term)
        } else {
          glm_model <- glm_model
          tested <- c(tested, max_pval$term)
          confounder <- c(confounder, max_pval$term)
        }

      }
    }
  }
  model_variables <- subset(glm_model_tbl$variable, glm_model_tbl$variable != "(Intercept)")
  results <- list(glm_model = glm_model, glm_model_tbl =  glm_model_tbl, model_variables = model_variables, confounder = confounder)
  return(results)
}

purposeful_selection_step_3 <-  function(outcome, fixed_model_variables, candidate_variables, dataset, retention_criteria, confounding_criteria) {
  # Input formula
  input_formula <- paste0(outcome, " ~ 1 + ", paste0(fixed_model_variables, collapse = "+"))
  input_formula <- trimws(input_formula, whitespace = "\\+")
  # Determine candidates for retension analysis
  model_additions <- lapply(candidate_variables, FUN = function(x) {
    glm_model <- glm(formula = paste0(input_formula, "+", x), data = dataset, family = "binomial")
    format_logistic_regression_table(glm_model)
 })
  # Curate model additions table
  model_additions_tbl <- do.call(rbind, model_additions)
  model_additions_tbl <- subset(model_additions_tbl, !variable %in% c("(Intercept)", fixed_model_variables))
  model_additions_tbl$p_value <- as.numeric(model_additions_tbl$p_value)
  model_additions_tbl <- arrange(model_additions_tbl, -p_value)
  # Final candidates
  final_candidates <- model_additions_tbl[model_additions_tbl[["p_value"]] < retention_criteria, "variable"]
  if (length(final_candidates) == 0) {
    glm_model <- glm(formula = input_formula, data = dataset, family = "binomial")
    glm_model_tbl <- format_logistic_regression_table(glm_model)
  } else {
    updated_formula <- paste0(input_formula, "+", paste0(final_candidates, collapse = "+"))
    glm_model <- glm(formula = updated_formula, data = dataset, family = "binomial")
    glm_model_tbl <- format_logistic_regression_table(glm_model)
    tested <- c()
    # Iterate through variables (i.e., test if p-value >1)
    repeat {
      # grab significant value
      if (length(tested) == length(final_candidates)) {
        break
      }

      eligible_variables_tbl <- glm_model_tbl[!glm_model_tbl[["variable"]] %in% c(tested, fixed_model_variables, "(Intercept)"), ]
      max_pval <- eligible_variables_tbl[which.max(eligible_variables_tbl[["p_value"]]), ]
      if (max_pval$p_value < retention_criteria) {
        glm_model <- glm_model
        tested <- c(tested, max_pval$term)
      } else {
        # Test as confounder
        ## Remove variable from model
        # Test confounding
        confounding <- test_confounding(variable = max_pval$variable, glm_model = glm_model, dataset = dataset, confounding_criteria = confounding_criteria)
        if (confounding$confounder_status == FALSE) {
          # Remove variable from model
          glm_model <- confounding[["reduced_model"]]
          # Regenerate summary
          glm_model_tbl <- format_logistic_regression_table(glm_model)
          # Label tested variables
          tested <- c(tested, max_pval$term)
        } else {
          glm_model_tbl <- format_logistic_regression_table(glm_model)
          tested <- c(tested, max_pval$variable)
        }
      }

    }
  }
  results <- list(candidate_variables = candidate_variables, model_additions_tbl = model_additions_tbl, final_model = glm_model, final_model_tbl = glm_model_tbl)
  return(results)
}

test_confounding <- function(variable, glm_model, dataset, confounding_criteria) {
  # Remove term from model
  reduced_model_formula <- gsub(variable, "", glm_model$formula)
  reduced_model_formula <- gsub("\\+\\+", "\\+", reduced_model_formula)
  reduced_model_formula <- trimws(reduced_model_formula, whitespace = "\\+")
  # Reduced model
  glm_reduced_model <- glm(formula=reduced_model_formula, data = dataset, family = "binomial")
  # Generate variable with
  full_model_summary <-  format_logistic_regression_table(glm_model)[,c("variable", "OR")]
  colnames(full_model_summary) <- c("variable", "OR_full")
  reduced_model_summary <- format_logistic_regression_table(glm_reduced_model)[,c("variable", "OR")]
  colnames(reduced_model_summary) <- c("variable", "OR_reduced")
  # Compare models
  model_comparison <- suppressMessages(left_join(full_model_summary, reduced_model_summary))
  model_comparison <- model_comparison[model_comparison[["variable"]] != tested_variable, ]
  # Convert OR to numeric
  model_comparison[["OR_full"]] <- as.numeric(model_comparison[["OR_full"]] )
  model_comparison[["OR_reduced"]] <- as.numeric(model_comparison[["OR_reduced"]] )
  # Test for confounding
  model_comparison$effect_change <- 1 - (model_comparison[["OR_full"]] / model_comparison[["OR_reduced"]])
  model_comparison$confounder <- ifelse(abs(model_comparison[["effect_change"]]) > confounding_criteria, "yes", "no")
  # IS it a confounder
  confounder_status <- ifelse(sum(model_comparison$confounder == "yes") > 0, TRUE, FALSE)
  # Results
  results <- list(confounder_status = confounder_status, confounder_tbl = model_comparison, reduced_model = glm_reduced_model)
  return(results)
}
