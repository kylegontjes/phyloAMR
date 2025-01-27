asr <- function(df,tr,tip_name_var ,pheno,model="ER",node_states = "marginal",conf_threshold=0.875){
  # Info

  library(tidyverse)
  library(corHMM)

   # Run corHMM to estimate hidden rates
  corHMM_out = corHMM(phy=tr,data=df[,c(tip_name_var,pheno)],rate.cat = 1,model=model,node.states = node_states)

  # Get Parent Child data
  outcome_str <- df[,pheno] %>% `names<-`(df[[tip_name_var]])
  parent_child_df <- get_parent_child_data(tr=tr,anc_data=corHMM_out$states,pheno_data=outcome_str,conf_threshold = conf_threshold)

  # Annotate parent child data
  parent_child_df <- get_phenotypic_continuation_data(parent_child_df)

  return(parent_child_df)
}

get_parent_child_data <- function(tr,anc_data,pheno_data,conf_threshold=0.875){
  # Tree edge info
  anc_data = as.data.frame(anc_data) %>% `colnames<-`(c("0","1"))
  edge <- tr$edge %>% as.data.frame
  colnames(edge) <- c("parent", "child")

  # Prediction data
  anc_data[,"node"] <-  1:nrow(anc_data)  + (nrow(anc_data)+1)
  anc_data[,"pred"] <- ifelse(anc_data[,2] > conf_threshold,1,ifelse(anc_data[,1] >conf_threshold,0,0.5))
  num_edges <- nrow(edge)

  # Dataframe construction
  edge$parent_mle <- rep(NA,num_edges)
  edge$parent_val <- rep(NA,num_edges)
  edge$child_mle <- rep(NA,num_edges)
  edge$child_val <- rep(NA,num_edges)

  # Prediction values for the parent & child
  for(i in nodes){
    edge[edge$parent == i,"parent_mle"] <- anc_data[anc_data$node == i,'1'] *100
    edge[edge$parent == i,"parent_val"] <- anc_data[anc_data$node == i,"pred"]
    edge[edge$child == i,"child_mle"] <- anc_data[anc_data$node == i,'1'] *100
    edge[edge$child == i,"child_val"] <- anc_data[anc_data$node == i,"pred"]
  }

  # Values and name for the tip
  for(i in 1:(min(nodes)-1)){
    edge[edge$child == i,"child_val"] <- pheno_data[[i]]
    edge[edge$child == i,"child_name"] <- names(pheno_data[i])
  }

  return(edge)
}

get_phenotypic_continuation_data <- function(parent_child_df){
  parent_child_df <- parent_child_df %>% dplyr::mutate(transition_any = ifelse(parent_val != child_val,1,0),
                                                       transition_high = ifelse(parent_val ==0 & child_val==1 | parent_val==1 & child_val ==0,1,0),
                                                       transition_low = ifelse(parent_val ==0.5 & child_val==1 | parent_val==0.5 & child_val ==0,1,0),
                                                       gain_low = ifelse(child_val ==1 & parent_val==0.5,1,0),
                                                       gain_high  =  ifelse(child_val ==1 & parent_val==0,1,0),
                                                       gain_any = ifelse(gain_low==1 | gain_high==1,1,0),
                                                       loss_low  =  ifelse(child_val ==0 & parent_val==0.5,1,0) ,
                                                       loss_high  =  ifelse(child_val ==0 & parent_val==1,1,0) ,
                                                       loss_any = ifelse(loss_low==1 | loss_high==1,1,0),
                                                       continuation_any  =  ifelse(parent_val == child_val,1,0),
                                                       continuation_high  =  ifelse(parent_val == 1 & child_val ==1 | parent_val == 0 & child_val ==0,1,0),
                                                       continuation_low =  ifelse(parent_val == 0.5 & child_val ==0.5,1,0))
  return(parent_child_df)
}
