find_best_model = function(df,tr,tip_name_var,pheno,node_states = "joint"){
  # Check if phenotype is 0,1
  check_phenotype(df[[pheno]])

  # Run corHMM to estimate hidden rates
  corHMM_ER = corHMM::corHMM(phy=tr,data=df[,c(tip_name_var,pheno)],rate.cat = 1,model="ER",node.states = node_states)
  corHMM_ARD = corHMM::corHMM(phy=tr,data=df[,c(tip_name_var,pheno)],rate.cat = 1,model="ARD",node.states = node_states)
  corHMM_SYS = corHMM::corHMM(phy=tr,data=df[,c(tip_name_var,pheno)],rate.cat = 1,model="SYS",node.states = node_states)

  AICS = c("ER"=corHMM_ER$AIC,"ARD"=corHMM_ARD$AIC,"SYS" = corHMM_SYS$AIC) %>% sort

  best_model = names(AICS)[corHMM_ER$AICc %>% which.min]
  return(best_model)

}
