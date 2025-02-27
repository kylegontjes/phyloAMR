#' Characterize binary ancestral state reconstruction model
#'
#' This function provides general model statistics for a binary ancestral state reconstruction model 
#'
#' @param corHMM_obj Ancestral reconstruction model from corHMM::corHMM 
#' @return A dataframe containing the log-likelihoods of the models, the rates and standard errors, and the model type  
#' @export
characterize_asr_model <- function(corHMM_obj){
  model <- ifelse(c(corHMM_obj[["index.mat"]] %>% as.vector() %>% subset(is.na(.)==F) %>% unique %>% length)==1,"ER","ARD")
  rate1 <- corHMM_obj$solution[1,2] %>% round(.,2)  
  rate2 <-  ifelse(model=="ARD",corHMM_obj$solution[2,1] %>% round(.,2)  ,NA)
  np = max(unlist(corHMM_obj$index.mat),na.rm=T)
  nRateCat = corHMM_obj$rate.cat
  loglik = corHMM_obj$loglik
  AIC = corHMM_obj$AIC
  AICc = corHMM_obj$AICc
  data <- data.frame(model,np,nRateCat,rate1,rate2,loglik,AIC,AICc)
  return(data)
}
