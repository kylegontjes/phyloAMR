#' phyloaware_regression: Phylogenetically-aware regression
#'
#' This function performs regression on our dataset
#'
#' @param trait Trait of interest
#' @param variables Exposure variables of interest
#' @param df Dataframe
#' @param first_present Whether to identify first present isolate from a patient
#' @param patient_id Patient identifier variable stored in dataframe df
#' @param culture_date Culture date variable stored in dataframe df
#' @param multivariable Boolean (i.e., TRUE, FALSE) whether to perform multivariable regression
#' @param stepwise_direction Direction if stepwise multivariable regression is chosen
#' @param entry_criteria P-value for defining candidate variables for multivariable regression (used in pvalue and purposeful selection)
#' @param retention_criteria P-value for retaining candidate variables in model (used in pvalue and purposeful selection)
#' @param confounding_criteria Impact on effect size for purposeful selection (used in purposeful selection)
#' @return List with univariable and multivariable results (if requested)
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats confint
#' @importFrom stats filter
#' @importFrom stats glm
#' @importFrom stats na.omit
#' @importFrom stats setNames
#' @export
phyloaware_regression <- function(trait, variables, df, first_present = NULL, patient_id = NULL, culture_date = NULL, multivariable = NULL, stepwise_direction = NULL, entry_criteria = NULL, retention_criteria = NULL, confounding_criteria = NULL) {
  # Get dataset with first isolate
  datasets <- phyloaware_dataset_curation(trait = trait, df = df, first_present = first_present, patient_id = patient_id, culture_date = culture_date)
  # Unadjusted
  univariable <- lapply(datasets, FUN = function(x) {
    univariable_regression_table(outcome = trait, dataset = x, variables = variables)
  })  %>% `names<-`(names(datasets))
  results <- list(datasets = datasets, univariable = univariable)
  # Multivariable
  if (is.null(multivariable) == TRUE | multivariable == FALSE) {
    return(results)
  } else {
    if (multivariable == "purposeful") {
    multivariable <- lapply(datasets, FUN = function(x) {
      purposeful_selection_algorithm(outcome = trait, variables = variables, dataset = x, entry_criteria = entry_criteria, retention_criteria = retention_criteria, confounding_criteria = confounding_criteria)
    })
  } else if (multivariable == "AIC") {
    multivariable <- lapply(datasets, FUN = function(x) {
      AIC_stepwise_regression(outcome = trait, dataset = x, variables = variables, stepwise_direction = stepwise_direction)
    })
  } else if (multivariable == "pvalue") {
    multivariable <- lapply(datasets, FUN = function(x) {
      pvalue_informed_regression(outcome = trait, dataset = x, variables = variables, entry_criteria = entry_criteria, retention_criteria = retention_criteria)
  })
  } else if (multivariable == "multivariable") {
    multivariable <- lapply(datasets, FUN = function(x) {
      model <- as.formula(paste0(trait, " ~ 1 +", paste0(variables, collapse = " + ")))
      glm_model <- glm(model, data = x, family = "binomial")
      ci <- suppressMessages(confint(glm_model))
      final <- cbind(exp(cbind(OR = coef(glm_model), ci)) %>% round(., 2) %>% formatC(., format = "f", digits = 2),
                     abs(summary(glm_model)$coefficients[, "Pr(>|z|)"]) %>%
                       round(., 4) %>% formatC(., format = "f", digits = 4)) %>% subset(rownames(.) != "(Intercept)") %>%
        `colnames<-`(c("OR", "2.5%", "97.5%", "p_value")) %>%
        as.data.frame %>% mutate(`OR (95% CI)` = paste0(OR, " (", `2.5%`, "-", `97.5%`, ")")) %>% select(`OR (95% CI)`, p_value)
      return(final)
    })
  }
    names(multivariable) <- names(datasets)
    results[["multivariable"]] <- multivariable
  }
  return(results)
}
