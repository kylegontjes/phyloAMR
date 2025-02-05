asr <- function(df,tr,tip_name_var ,pheno,model="ER",node_states = "joint",conf_threshold=0.875){
  # Check if phenotype is 0,1
  check_phenotype(df[[pheno]])

  # Order dataframe properly
  df<- df %>% .[match(tr$tip.label,.[[tip_name_var]]),]

  # Check if you want to use the best model
  if(model=="MF"){
    print("Performing model finding to use best fitting model")
    model = find_best_model(df=df,tr=tr,tip_name_var=tip_name_var,pheno=pheno,node_states = node_states)
    print(paste0("Best model is: ",model))
  }

  # Run corHMM to estimate hidden rates
  corHMM_out = corHMM::corHMM(phy=tr,data=df[,c(tip_name_var,pheno)],rate.cat = 1,model=model,node.states = node_states)

  # Get Parent Child data
  outcome_str <- df[,pheno] %>% `names<-`(df[[tip_name_var]])
  parent_child_df <- get_parent_child_data(tr=tr,anc_data=corHMM_out$states,pheno_data=outcome_str,conf_threshold = conf_threshold,node_states=node_states)

  # Annotate parent child data
  parent_child_df <- get_phenotypic_continuation_data(parent_child_df)

  asr_output = list(corHMM_out = corHMM_out,parent_child_df = parent_child_df)
  return(asr_output)
}

check_phenotype <- function(phenotype_var){
  if(sum(unique(as.numeric(phenotype_var)) %in% c(0,1))!=2){
    stop("Phenotype is not formatted as binary variable with event = 1 and no-event = 0")
  }
}

get_parent_child_data <- function(tr,anc_data,pheno_data,conf_threshold=0.875,node_states){
  # Tree edge info
  edge <- tr$edge %>% as.data.frame
  colnames(edge) <- c("parent", "child")

  # Prediction data
  if(node_states == "marginal"){
    anc_data = as.data.frame(anc_data) %>% `colnames<-`(levels(as.factor(pheno_data)))
    anc_data[,"pred"] <- ifelse(anc_data[,2] > conf_threshold,1,ifelse(anc_data[,1] >conf_threshold,0,0.5))
  } else if (node_states == "joint"){
    state_coding = as.factor(pheno_data) %>% levels %>% `names<-`(1:length(.))
    anc_data = data.frame(pred = anc_data)
    anc_data$pred =  recode(anc_data$pred,!!!as.list(state_coding)) %>% as.numeric
  } else {
    stop("Phylosuite is only optimized for marginal and joint corHMM")
  }

  # Add node data
  anc_data[,"node"] <-  1:nrow(anc_data)  + (nrow(anc_data)+1)
  num_edges <- nrow(edge)
  internal_nodes <- c(length(pheno_data) +1):(unlist(edge) %>% max)

  # Dataframe construction
  edge$parent_val <- rep(NA,num_edges)
  edge$child_val <- rep(NA,num_edges)

  # Prediction values for the parent & child
  for(i in internal_nodes){
    edge[edge$parent == i,"parent_val"] <- anc_data[anc_data$node == i,"pred"]
    edge[edge$child == i,"child_val"] <- anc_data[anc_data$node == i,"pred"]
  }

  # Internal nodes
  if(node_states == "marginal"){
    edge$child_mle <- rep(NA,num_edges)
    edge$parent_mle <- rep(NA,num_edges)
    for(i in internal_nodes){
    edge[edge$parent == i,"parent_mle"] <- anc_data[anc_data$node == i,levels(as.factor(pheno_data))[2]] *100
    edge[edge$child == i,"child_mle"] <- anc_data[anc_data$node == i,levels(as.factor(pheno_data))[2]] *100
    }
  }

  # Values and name for the tip
  for(i in 1:(min(internal_nodes)-1)){
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

