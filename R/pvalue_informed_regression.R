pvalue_informed_regression <- function(outcome,dataset,variables,entry_criteria,retention_criteria){
  datatable <- lapply(variables,FUN=function(x){
    datatable <- univariable_regression(x,outcome,dataset) %>% as.data.frame
    return(datatable)
  })%>% do.call(rbind,.)

  pval_less_threshold <- datatable %>% subset(p_value < as.numeric(entry_criteria))  %>% arrange(p_value)
  variables <- rownames(pval_less_threshold)
  # Loop through variables to create final model
  if(min(pval_less_threshold$p_value)>retention_criteria){
    return(paste0("No model - no variables < ",retention_criteria))
  }
  if(min(pval_less_threshold$p_value)<retention_criteria){
    # Now loop through variables
    null_model <- paste(outcome,"1", sep = "~")
    null_model_out <- glm(null_model, data = dataset, family = "binomial")
    # Retain if pvalue is
    best_variable <- variables[1]
    null_model <- paste0(null_model_out$formula,"+",best_variable)
    null_model_out <- glm(null_model, data = dataset, family = "binomial")
    for(i in variables[2:length(variables)]){
      null_model_new <- paste0(null_model_out$formula,"+",i)
      null_model_out_new <- glm(null_model_new, data = dataset, family = "binomial")
      pvals <- abs(summary(null_model_out_new)$coefficients[,'Pr(>|z|)']) %>% round(.,4)
      if(pvals[length(pvals)] >retention_criteria){
        null_model_out <- null_model_out
      }
      if(pvals[length(pvals)] <retention_criteria){
        null_model_out <- null_model_out_new
      }
    }
    output <- null_model_out
    ci <- suppressMessages(confint(output))
    final <- cbind(exp(cbind(OR = coef(output), ci)) %>% round(.,2),abs(summary(output)$coefficients[,'Pr(>|z|)']) %>% round(.,4)) %>% subset(rownames(.) != "(Intercept)")  %>% `colnames<-`(c("OR","2.5%","97.5%","p_value")) %>%  as.data.frame %>%
      mutate(`OR (95% CI)` = paste0(OR," (",`2.5%`,"-",`97.5%`,")")) %>% select(`OR (95% CI)`,p_value)
    results <- list(final_model=final,univariable=datatable)
    return(results)
  }
}
