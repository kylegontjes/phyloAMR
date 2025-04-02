#' phyloaware_regression: Phylogenetically-aware regression
#'
#' This function performs regression on our dataset
#'
#' @param pheno Trait of interest
#' @param variables Exposure variables of interest
#' @param df Dataframe
#' @param first_present Whether to identify first present isolate from a patient
#' @param patient_id Patient identifier variable stored in dataframe df
#' @param culture_date Culture date variable
#' @param multivariable Boolean whether to perform multivariable regression
#' @param stepwise_direction Direction if stepwise multivariable regression is chosen
#' @param entry_criteria P-value for defining candidate variables for multivariable regression
#' @param retention_criteria P-value for retaining candidate variables in model
#' @param confounding_criteria Impact on effect size
#' @return List with univariable and multivariable results (if requested)
#' @export
phyloaware_regression <- function(pheno, variables, df, first_present = NULL, patient_id = NULL, culture_date = NULL, multivariable = NULL, stepwise_direction = NULL, entry_criteria = NULL, retention_criteria = NULL, confounding_criteria = NULL) {
  # Get dataset with first isolate
  datasets <- phyloaware_dataset_curation(pheno, df, first_present = first_present, patient_id = patient_id, culture_date = culture_date)
  # Unadjusted
  univariable <- lapply(datasets, FUN = function(x) {
    univariable_regression_table(outcome = pheno, dataset = x, variables = variables)
  })  %>% `names<-`(names(datasets))
  results <- list(datasets = datasets, univariable = univariable)
  # Multivariable
  if (is.null(multivariable) == TRUE) {
    stop("No multivariable regression requested")
    return(results)
  }
  if (multivariable == "purposeful") {
    multivariable <- lapply(datasets, FUN = function(x) {
      purposeful_selection_algorithm(outcome = pheno, variables = variables, dataset = x, entry_criteria = entry_criteria, retention_criteria = retention_criteria, confounding_criteria = confounding_criteria)
    })  %>% `names<-`(names(datasets))
    results[["multivariable"]] <- multivariable
  } else if (multivariable == "AIC") {
    multivariable <- lapply(datasets, FUN = function(x) {
      model <- as.formula(paste0(pheno, " ~ 1 +", paste0(variables, collapse = " + ")))
      if (stepwise_direction == "forward") {
        initialize_model <- glm(formula =  as.formula(paste0(pheno, " ~ 1")), data = x, family = "binomial")
        stepwise_glm_model <- suppressMessages(stats::step(initialize_model, direction = stepwise_direction, scope = model))
      }   else {
        glm_model <- glm(model, data = x, family = "binomial")
        stepwise_glm_model <- suppressMessages(stats::step(glm_model, direction = stepwise_direction))
      }
      ci <- suppressMessages(confint(stepwise_glm_model))
      final <- cbind(exp(cbind(OR = coef(stepwise_glm_model), ci)) %>% round(., 2),
                     abs(summary(stepwise_glm_model)$coefficients[, "Pr(>|z|)"]) %>%
                       round(., 4)) %>% subset(rownames(.) != "(Intercept)") %>%
        `colnames<-`(c("OR", "2.5%", "97.5%", "p_value")) %>%
        as.data.frame %>% mutate(`OR (95% CI)` = paste0(OR, " (", `2.5%`, "-", `97.5%`, ")")) %>% select(`OR (95% CI)`, p_value)
      return(final)
    }) %>% `names<-`(names(datasets))
    results[["multivariable"]] <- multivariable
  } else if (multivariable == "pvalue") {
    multivariable <- lapply(datasets, FUN = function(x) {
      pvalue_informed_regression(outcome = pheno, dataset = x, variables = variables, pvalue_threshold = entry_criteria) %>% .["final_model"]
  }) %>% `names<-`(names(datasets))
    results[["multivariable"]] <- multivariable
  }
  return(results)
}
