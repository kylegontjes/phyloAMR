purposeful_selection_algorithm <- function(outcome,variables,dataset,entry_criteria,retention_criteria,confounding_criteria){
  ps_step1 <- purposeful_selection_step_1(outcome = outcome,variables = variables,dataset = dataset,entry_criteria = entry_criteria)
  ps_step2 <- purposeful_selection_step_2(outcome = outcome,candidate_variables = ps_step1$candidates,dataset = dataset,retention_criteria = retention_criteria,confounding_criteria = confounding_criteria)
  if(length(ps_step2)== 1 & sum(ps_step2 == "No significant variables")>0){
    results <- list(ps_step1 = ps_step1,warning = 'No multivariable model possible as all candidate variables were not statistically significant.')
  } else {
    p3_candidate_variables <- subset(variables,!variables %in% ps_step1$candidates)
    ps_step3 <- purposeful_selection_step_3(outcome = outcome,fixed_model_variables= ps_step2$model_variables,candidate_variables = p3_candidate_variables,dataset = dataset,retention_criteria = retention_criteria,confounding_criteria=confounding_criteria)
    final_model_table <- purposeful_table_curation(ps_step3$final_model)
    results <- list(ps_step1=ps_step1,ps_step2=ps_step2,ps_step3=ps_step3,final_model=ps_step3$final_model,final_model_table=final_model_table)
  }
  return(results)
}

purposeful_selection_step_1 <- function(outcome,variables,dataset,entry_criteria){
  univariable_results <- lapply(variables,FUN=function(x){glm(formula = paste0(outcome," ~ 1 +",x),data=dataset,family = "binomial")})
  univariable_results_tbl <- univariable_results %>% lapply(.,FUN=function(x){data.table::data.table(coef(summary(x)),keep.rownames='term')}) %>% do.call(rbind,.) %>% subset(term != '(Intercept)') %>% arrange(`Pr(>|z|)`)
  candidates <- subset(univariable_results_tbl,`Pr(>|z|)`<entry_criteria) %>% .$term
  non_candidates <- subset(variables,!variables %in% candidates)
  results <- list(candidates=candidates,non_candidates=non_candidates,univariable_results_tbl=univariable_results_tbl)
  return(results)
}

purposeful_selection_step_2 <- function(outcome,candidate_variables,dataset,retention_criteria,confounding_criteria){
  input_formula <- paste0(outcome," ~ 1 + ",paste0(candidate_variables,collapse="+")) %>% trimws(.,whitespace = "\\+")
  # Glm model output to manipulate
  glm_model <- glm(formula=input_formula,data=dataset,family="binomial")
  glm_model_tbl <- data.table::data.table(coef(summary(glm_model)),keep.rownames='term')
  # Variables tested
  tested <- c()
  confounder <- c()
  # Iterative process: For top model component, drop if alpha lvel is above 0.1 AND not a confounder (two checks)
  repeat{
    # grab significant value
    max_pval <- glm_model_tbl %>% subset(.,!term %in% c(tested,"(Intercept)"))  %>% .[which.max(.$`Pr(>|z|)`),]

    if(max_pval$`Pr(>|z|)` < retention_criteria){
      # Model is completed as all either confounders or significant
      break
    } else {
      if(c(length(tested)+1)==length(candidate_variables) & max_pval$`Pr(>|z|)` > retention_criteria){
        warning("No significant variables, reporting step 1 analysis")
        results <- "No significant variables"
        break
      } else {
        # Test as confounder
        ## Remove variable from model
        reduced_model_formula <- glm_model$formula %>% gsub(max_pval$term,"",.)  %>% gsub("\\+\\+","\\+",.) %>% trimws(.,whitespace = "\\+")
        glm_reduced_model <- glm(formula=reduced_model_formula,data = dataset,family="binomial")
        confounding <- test_confounding(model1=glm_model,model2=glm_reduced_model,confounding_criteria=confounding_criteria,tested_variable = max_pval$term)
        is_confounding <- ifelse(sum(confounding$confounder=="yes")>0,"yes","no")
        if(is_confounding =="no"){
          # Remove variable from model
          glm_model <- glm_reduced_model
          # Regenerate summary
          glm_model_tbl <- data.table::data.table(coef(summary(glm_model)),keep.rownames='term')
          # Label tested variables
          tested <- c(tested,max_pval$term)
        } else{
          glm_model <- glm_model
          tested <- c(tested,max_pval$term)
          confounder <- c(confounder,max_pval$term)

        }
        model_variables <- glm_model_tbl$term %>% subset(.!="(Intercept)")
        results <- list(glm_model=glm_model,glm_model_tbl=glm_model_tbl,model_variables =model_variables,confounder=confounder)
      }
    }
  }

  return(results)
}

purposeful_selection_step_3 <-  function(outcome,fixed_model_variables,candidate_variables,dataset,retention_criteria,confounding_criteria){
  input_formula <- paste0(outcome," ~ 1 + ",paste0(fixed_model_variables,collapse="+")) %>% trimws(.,whitespace = "\\+")
  # Determine candidates for retension analysis
  model_additions <- lapply(candidate_variables,FUN=function(x){glm(formula = paste0(input_formula,"+",x),data=dataset,family = "binomial")})
  model_additions_tbl <- model_additions %>% lapply(.,FUN=function(x){data.table::data.table(coef(summary(x)),keep.rownames='term')}) %>% do.call(rbind,.) %>% arrange(-`Pr(>|z|)`) %>% subset(!term  %in%c('(Intercept)',fixed_model_variables))
  final_candidates <- subset(model_additions_tbl,`Pr(>|z|)`<retention_criteria)   %>% .$term
  if(length(final_candidates)==0){
    glm_model <- glm(formula=input_formula,data=dataset,family="binomial")
    glm_model_tbl <- data.table::data.table(coef(summary(glm_model)),keep.rownames='term')
  } else {
    updated_formula <- paste0(input_formula,"+",paste0(final_candidates,collapse="+"))
    glm_model <- glm(formula=updated_formula,data=dataset,family="binomial")
    glm_model_tbl <- data.table::data.table(coef(summary(glm_model)),keep.rownames='term')
    tested <- c()
    # Iterate through variables (i.e., test if p-value >1)
    repeat{
      # grab significant value
      if(length(tested)==length(final_candidates)){
        break
      }
      max_pval <- glm_model_tbl %>% subset(.,!term %in% c(tested,fixed_model_variables,"(Intercept)"))  %>% .[which.max(.$`Pr(>|z|)`),]
      if(max_pval$`Pr(>|z|)` < retention_criteria){
        glm_model <- glm_model
        tested <- c(tested,max_pval$term)
      } else {
        # Test as confounder
        ## Remove variable from model
        reduced_model_formula <- glm_model$formula %>% gsub(max_pval$term,"",.)  %>% gsub("\\+\\+","\\+",.) %>% trimws(.,whitespace = "\\+")
        glm_reduced_model <- glm(formula=reduced_model_formula,data = dataset,family="binomial")
        confounding <- test_confounding(model1=glm_model,model2=glm_reduced_model,confounding_criteria=confounding_criteria,tested_variable = max_pval$term)
        is_confounding <- ifelse(sum(confounding$confounder=="yes")>0,"yes","no")
        if(is_confounding =="no"){
          # Remove variable from model
          glm_model <- glm_reduced_model
          # Regenerate summary
          glm_model_tbl <- data.table::data.table(coef(summary(glm_model)),keep.rownames='term')
          # Label tested variables
          tested <- c(tested,max_pval$term)
        } else{
          glm_model <- glm_model
          glm_model_tbl <- data.table::data.table(coef(summary(glm_model)),keep.rownames='term')
          tested <- c(tested,max_pval$term)

        }
      }

    }
  }
  results <- list(candidate_variables=candidate_variables,model_additions_tbl=model_additions_tbl,final_model=glm_model,final_model_tbl=glm_model_tbl)
  return(results)
}

test_confounding <- function(model1,model2,confounding_criteria,tested_variable){
  model1_summary <-  data.table::data.table(coef(summary(model1)),keep.rownames='term') %>% `colnames<-`(c("term","m1_estimate","m1_SE","m1_z_val","m1_p_val"))
  model2_summary <- data.table::data.table(coef(summary(model2)),keep.rownames='term')%>% `colnames<-`(c("term","m2_estimate","m2_SE","m2_z_val","m2_p_val"))
  model_comparison <- suppressMessages(left_join(model1_summary,model2_summary))
  model_comparison <- model_comparison  %>% subset(term !=tested_variable)
  model_comparison$effect_change <- 1 - (model_comparison$m1_estimate / model_comparison$m2_estimate)
  model_comparison$confounder <- ifelse(abs(model_comparison$effect_change)>confounding_criteria,"yes","no")
  return(model_comparison)
}

purposeful_table_curation <- function(final_model){
  ci <- suppressMessages(confint(final_model))
  final <- cbind(exp(cbind(OR = coef(final_model), ci)) %>% round(.,
                                                                  2), abs(summary(final_model)$coefficients[, "Pr(>|z|)"]) %>%
                   round(., 4)) %>% subset(rownames(.) != "(Intercept)") %>%
    `colnames<-`(c("OR", "2.5%", "97.5%", "p_value")) %>%
    as.data.frame %>% mutate(`OR (95% CI)` = paste0(OR, " (",
                                                    `2.5%`, "-", `97.5%`, ")")) %>% select(`OR (95% CI)`,
                                                                                           p_value)
  final <- final %>% select(`OR (95% CI)`,`p_value`)
  return(final)
}
