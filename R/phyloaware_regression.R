phyloaware_regression <- function(df,pheno,asr_cluster_obj,datasets,variables,entry_criteria,retention_criteria,confounding_criteria){
  # Get dataset with first isolate
  datasets= dataset_curation(pheno,df)
  # Purposeful selection
  lr_non_s <- tryCatch(purposeful_selection_algorithm(outcome = phenotype_var,dataset = datasets[[paste0(drug,"_non_s")]],variables = variables,entry_criteria = entry_criteria,retension_criteria=retension_criteria,confounding_criteria=confounding_criteria), error = function(e) e)
  lr_emerge <-  tryCatch(purposeful_selection_algorithm(outcome = phenotype_var,dataset = datasets[[paste0(drug,"_emerge")]],variables = variables,entry_criteria = entry_criteria,retension_criteria=retension_criteria,confounding_criteria=confounding_criteria), error = function(e) e)
  lr_spread <- tryCatch(purposeful_selection_algorithm(outcome = phenotype_var,dataset = datasets[[paste0(drug,"_spread")]],variables = variables,entry_criteria = entry_criteria,retension_criteria=retension_criteria,confounding_criteria=confounding_criteria), error = function(e) e)
  lr_list <- list(lr_non_s,lr_emerge,lr_spread) %>% `names<-`(c("non_s","emerge","spread"))
  return(lr_list)
}
