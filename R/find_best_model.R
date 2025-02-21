#' find_best_model: Determine best model for ancestral state reconstruction
#'
#' Function to find best model for ancestral state reconstruction using corHMM.
#'
#' @param df Dataframe with tip name variable and phenotype
#' @param tr Phylogenetic tree
#' @param tip_name_var Name of variable containing tip names in df
#' @param pheno Name of phenotype variable in df
#' @param node_states Whether to perform "joint" or "marginal" reconstruction
#' @return Description of return value
#'   \describe{
#'     \item{best_model}{Best model to use for the ancestral state reconstruction}
#'   }
#' @export

find_best_model = function(df,tr,tip_name_var,pheno,node_states = "joint",upper_bound=1e9,lower_bound=1e-9){
  # Check if phenotype is 0,1
  check_phenotype(df[[pheno]])

  # Run corHMM to estimate hidden rates
  corHMM_ER = corHMM::corHMM(phy=tr,data=df[,c(tip_name_var,pheno)],rate.cat = 1,model="ER",node.states = node_states,upper.bound=upper_bound,lower.bound=lower_bound)
  corHMM_ARD = corHMM::corHMM(phy=tr,data=df[,c(tip_name_var,pheno)],rate.cat = 1,model="ARD",node.states = node_states,upper.bound=upper_bound,lower.bound=lower_bound)
  corHMM_SYS = corHMM::corHMM(phy=tr,data=df[,c(tip_name_var,pheno)],rate.cat = 1,model="SYS",node.states = node_states,upper.bound=upper_bound,lower.bound=lower_bound)

  AICS = c("ER"=corHMM_ER$AIC,"ARD"=corHMM_ARD$AIC,"SYS" = corHMM_SYS$AIC) %>% sort

  best_model = names(AICS)[corHMM_ER$AICc %>% which.min]
  return(best_model)

}
