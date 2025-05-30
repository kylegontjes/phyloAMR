#' Dataset curation for phylogenetically aware regression
#'
#' Curate three datasetes for phylogenetically aware regression.
#'
#' Datasets included:
#' 1. Present: All isolates
#' 2. Singleton: Singleton isolates-only
#' 3. Cluster: Cluster isolates-only
#'
#' This permits isolated analysis of singleton and cluster isolates, enabling identification of characteristics associated with the emergence and/or spread of a trait.
#'
#' @param trait Outcome of interest. Character string.
#' @param df Dataset with trait variable and asr_cluster. Must contain patient_id and culture_date variables if first_present == TRUE.
#' @param first_present Boolean (i.e., TRUE/FALSE) indicating whether to take the first present isolate
#' @param patient_id Variable that contains IDs that group tips (i.e., Patient ID). Character string.
#' @param culture_date  Variable that contains collection/culture date. Character string. Variable must be formatted as a date.
#' @return Three datasets labeled as 'present', 'singleton', and 'cluster'.
#' @importFrom dplyr slice
#' @importFrom dplyr group_by
#' @importFrom dplyr across
#' @importFrom dplyr all_of
#' @importFrom dplyr arrange
#' @export

phyloaware_dataset_curation <- function(trait, df, first_present = FALSE, patient_id = NULL, culture_date = NULL) {
  # Check dataset has trait and asr
  check_asr_trait(df = df, trait = trait)

  # Get first isolate present
  if (first_present == TRUE) {
    # Check if patient_id and culture_date are present
    check_patient_id_and_culture_date(df = df, patient_id = patient_id, culture_date = culture_date)
    # Generate dataset with first present isolate
    present_df <- get_dataset_with_first_present_isolate(variable = trait, patient_id = patient_id, culture_date = culture_date, df = df)
  } else {
    present_df <- as.data.frame(df)
  }
  # Singleton-only dataset
  singleton_df <- subset(present_df, get(trait) == 0 | grepl("singleton|1pt", asr_cluster))
  # Cluster-only dataset
  cluster_df <-  subset(present_df, get(trait) == 0 | (grepl("cluster", asr_cluster) & !grepl("1pt", asr_cluster)))
  # Final list with overall, singleton, and cluster datasets
  df_list <- list(present_df, singleton_df, cluster_df)
  names(df_list) <- c("present", "singleton", "cluster")
  return(df_list)
}

get_dataset_with_first_present_isolate <- function(variable, patient_id, culture_date, df) {
  one_isolate <-  names(which(table(df[[patient_id]]) == 1))
  # Only 1 isolate
  df_one_isolate <- df[df[[patient_id]] %in% one_isolate, ]
  # Patient w/ more than one
  df_multiple <- df[!df[[patient_id]] %in% one_isolate, ]
  ## Determine if multiple resistant
  at_least_one_present <- unique(df_multiple[df_multiple[[variable]] == 1, patient_id])
  ## For only susceptible
  df_multiple_absent_only_first <- df_multiple[!df_multiple[[patient_id]] %in% at_least_one_present, ] %>% group_by(across(all_of(patient_id))) %>% arrange(across(all_of(culture_date))) %>% slice(1)
  ## Patients w/ >=1 resistant isolate
  df_multiple_first_resistant <-  df_multiple[df_multiple[[patient_id]] %in% at_least_one_present, ] %>% subset(get(variable) == 1) %>% group_by(across(all_of(patient_id))) %>% arrange(across(all_of(culture_date))) %>% slice(1)
  # Bind one isolate, susceptible first isolate, and first resistant isolate for paients w/ >=1 isolate
  first_present_isolate_df <- rbind(rbind(df_one_isolate, df_multiple_absent_only_first), df_multiple_first_resistant)
  return(first_present_isolate_df)
}
