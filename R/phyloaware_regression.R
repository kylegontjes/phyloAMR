#' phyloaware_regression: Phylogenetically-aware regression
#'
#' This function performs regression on our dataset
#'
#' @param trait Trait of interest
#' @param variables Exposure variables of interest
#' @param df Dataframe
#' @param first_present Boolean (i.e., TRUE, FALSE) whether to select a participant's first isolate with the trait.
#' @param patient_id Patient identifier variable stored in dataframe df. Required if first_present == TRUE.
#' @param culture_date Culture date used to select the participant's first isolate. Must be stored as a variable in the dataframe df. Required if first_present == TRUE.
#' @param multivariable Defines the multivariable selection strategy. Options include: 'purposeful', 'AIC', 'pvalue', and 'multivariable.'
#' @param stepwise_direction Direction of stepwise selection. Options include: 'both', 'backward', or 'forward'.
#' @param entry_criteria P-value for defining candidate variables for multivariable regression. Used in pvalue and purposeful selection)
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
  # Curate the three datasets (i.e., present, singletons, clusters)
  datasets <- phyloaware_dataset_curation(trait = trait, df = df, first_present = first_present, patient_id = patient_id, culture_date = culture_date)
  # Univariable regression
  univariable <- lapply(datasets, FUN = function(x) {
    univariable_regression(outcome = trait, dataset = x, variables = variables)
  })
  names(univariable) <- names(datasets)
  results <- list(datasets = datasets, univariable = univariable)
  # Multivariable
  if (is.null(multivariable) == TRUE | multivariable == FALSE) {
    return(results)
  } else {
    if (multivariable == "purposeful") {
    # Purposeful selection
    multivariable <- lapply(datasets, FUN = function(x) {
      purposeful_selection_algorithm(outcome = trait, variables = variables, dataset = x, entry_criteria = entry_criteria, retention_criteria = retention_criteria, confounding_criteria = confounding_criteria)
    })
  } else if (multivariable == "AIC") {
    # Stepwise regression using AIC
    multivariable <- lapply(datasets, FUN = function(x) {
      AIC_stepwise_regression(outcome = trait, dataset = x, variables = variables, stepwise_direction = stepwise_direction)
    })
  } else if (multivariable == "pvalue") {
    # P-value informed
    multivariable <- lapply(datasets, FUN = function(x) {
      pvalue_informed_regression(outcome = trait, dataset = x, variables = variables, entry_criteria = entry_criteria, retention_criteria = retention_criteria)
  })
  } else if (multivariable == "multivariable") {
    # Multivariable regression
    multivariable <- lapply(datasets, FUN = function(x) {
      multivariable_regression(outcome = trait, variables = variables, dataset = x)
      return(final)
    })
  }
    names(multivariable) <- names(datasets)
    results[["multivariable"]] <- multivariable
  }
  return(results)
}
