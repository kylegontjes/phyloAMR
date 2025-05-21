phyloaware_dataset_curation <- function(trait,df,first_present=NULL,patient_id=NULL,culture_date=NULL){
  # Get first isolate present
  if(first_present==T){
    present_df <- get_dataset_with_first_present_isolate(variable=trait, patient_id=patient_id, culture_date=culture_date,df=df)
  } else {
    present_df <- df %>% as.data.frame
  }
  # Emergence
  singleton_df <- present_df %>% subset(get(trait) ==0 | grepl("singleton|1pt",asr_cluster)) %>% as.data.frame()
  # Spread
  cluster_df <- present_df %>% subset(get(trait) ==0 | (grepl("cluster",asr_cluster) & !grepl("1pt",asr_cluster))) %>% as.data.frame()
  # Final DF list
  df_list <- list(present_df,singleton_df,cluster_df) %>% `names<-`(c("present","singleton","cluster"))
  return(df_list)
}

get_dataset_with_first_present_isolate <- function(variable,patient_id,culture_date,df){
  one_isolate <- df %>% group_by(across(patient_id)) %>% filter(n()==1) %>% .[[patient_id]]
  # Only 1 isolate
  df_one_isolate <- df[df[[patient_id]]%in% one_isolate,]
  # Patient w/ more than one
  df_multiple <- df[!df[[patient_id]]%in% one_isolate,]
  ## Determine if multiple resistant
  at_least_one_present <- unique(df_multiple[df_multiple[[variable]]==1,patient_id])
  ## For only susceptible
  absent_only_first <- df_multiple[!df_multiple[[patient_id]] %in% at_least_one_present,] %>% group_by(across(patient_id)) %>% arrange(get(culture_date)) %>% filter(row_number()==1)  %>% ungroup()
  ## Patients w/ >=1 resistant isolate
  first_resistant <-  df_multiple[df_multiple[[patient_id]] %in% at_least_one_present,] %>% group_by(across(patient_id)) %>% arrange(get(culture_date)) %>% subset(get(variable)==1) %>% filter(row_number()==1)  %>% ungroup()
  # Bind one isolate, susceptible first isolate, and first resistant isolate for paients w/ >=1 isolate
  first_present_isolate_df <- rbind(df_one_isolate,absent_only_first) %>% rbind(.,first_resistant)
  return(first_present_isolate_df)
}
