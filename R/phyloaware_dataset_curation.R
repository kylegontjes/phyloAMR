phyloaware_dataset_curation <- function(trait, df, first_present = NULL, patient_id = NULL, culture_date = NULL) {
  # Get first isolate present
  if (first_present == TRUE) {
    present_df <- get_dataset_with_first_present_isolate(variable = trait, patient_id = patient_id, culture_date = culture_date, df = df)
  } else {
    present_df <- as.data.frame(df)
  }
  # Emergence
  singleton_df <- subset(present_df, get(trait) == 0 | grepl("singleton|1pt", asr_cluster))
  # Spread
  cluster_df <-  subset(present_df, get(trait) == 0 | (grepl("cluster", asr_cluster) & !grepl("1pt", asr_cluster)))
  # Final DF list
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
