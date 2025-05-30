#' Phylogenetically aware regression of genome-influenced traits
#'
#' This function performs regression on our dataset
#'
#' Alongside univariable regression, multivariable options included in this implementation:
#' 1. Multivariable: Standard multivariable regression of all variables
#' 2. pvalue: A p-value informed logistic regression
#' 3. AIC: Stepwise AIC
#' 4. Purposeful selection: Iterative model selection accounting for p-value and confounding. Hosmer & Lemeshow 2000
#'
#' @param trait Trait of interest
#' @param variables Exposure variables of interest. Must be numeric or one-hot encoded
#' @param df Dataframe that contains the trait, exposure variables, and other requested variables.
#' @param first_present Boolean (i.e., TRUE, FALSE) whether to select a participant's first isolate with the trait.
#' @param patient_id Patient identifier variable stored in dataframe df. Required if first_present == TRUE.
#' @param culture_date Culture date used to select the participant's first isolate. Must be stored as a variable in the dataframe df. Required if first_present == TRUE.
#' @param multivariable Defines the multivariable selection strategy. Options include: 'purposeful', 'AIC', 'pvalue', and 'multivariable.'
#' @param stepwise_direction Direction of stepwise selection. Options include: 'both', 'backward', or 'forward'.  For more information, see stats::step.
#' @param entry_criteria P-value for defining candidate variables for multivariable regression. Used in pvalue and purposeful selection. Suggestion: 0.2.
#' @param retention_criteria P-value for retaining candidate variables in model Used in pvalue and purposeful selection. Suggestion: 0.1.
#' @param confounding_criteria Percent change of effect size. Set value to be very high (i.e., 1000) if testing for confounding is not desired. Default = 0.2.
#' @return List with datasets, univariable, and multivariable results (if requested)
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats confint
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
  results <- list(datasets = datasets, univariable = univariable)

  # Multivariable modeling
  if (is.null(multivariable) == TRUE || multivariable == FALSE) {
    return(results)
  } else {
    if (multivariable == "purposeful") {
      # Purposeful selection
      multivariable_results <- lapply(datasets, FUN = function(x) {
        purposeful_selection_algorithm(outcome = trait, variables = variables, dataset = x, entry_criteria = entry_criteria, retention_criteria = retention_criteria, confounding_criteria = confounding_criteria)
      })
    } else if (multivariable == "AIC") {
      # Stepwise regression using Akaike information criterion (AIC) via stats::step
      multivariable_results <- lapply(datasets, FUN = function(x) {
        AIC_stepwise_regression(outcome = trait, dataset = x, variables = variables, stepwise_direction = stepwise_direction)
      })
    } else if (multivariable == "pvalue") {
      # P-value informed regression
      multivariable_results <- lapply(datasets, FUN = function(x) {
        pvalue_informed_regression(outcome = trait, dataset = x, variables = variables, entry_criteria = entry_criteria, retention_criteria = retention_criteria)
      })
    } else if (multivariable == "multivariable") {
      # Multivariable regression without variable selection
      multivariable_results <- lapply(datasets, FUN = function(x) {
        multivariable_regression(outcome = trait, variables = variables, dataset = x)
      })
    }
    results[["multivariable"]] <- multivariable_results
  }
  return(results)
}
