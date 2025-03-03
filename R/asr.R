#' asr: ancestral state reconstruction wrapper function
#'
#' Function to perform ancestral state reconstruction using corHMM. Parse output into informative dataset with episodes of trait gain, loss, and continuation across the phylogenetic tree
#'
#' @param df Dataframe with tip name variable and phenotype
#' @param tr Phylogenetic tree
#' @param tip_name_var Name of variable containing tip names in df
#' @param pheno Name of phenotype variable in df
#' @param model Whether to use equal rates "ER" or all-rates differ "ARD" rate matrices.
#' @param node_states Whether to perform "joint" or "marginal" reconstruction
#' @param upper.bound Upper bound for likelihood search. The default for my implementation is 1e9
#' @param lower.bound Lower bound for likelihood search. The default for my implementation is 1e-9
#' @param conf_threshold The confidence threshold to use for marginal state reconstruction
#' @return Description of return value
#'   \describe{
#'     \item{corHMM_out}{corHMM output}
#'     \item{parent_child_df}{A dataframe of parent-child relationships with additional descriptive information.}
#'   }
#' @export
asr <- function(df,tr,tip_name_var ,pheno, model="ER", node_states = "joint", upper_bound=1e100, lower_bound=1e-9, conf_threshold=NULL){
  # Check if phenotype is 0,1
  check_phenotype(df[[pheno]])
  # Check if tree is rooted
  check_tree(tr)
  # Check if phenotype and tree length is the same
  check_phenotype_tree_length(pheno_var =df[[pheno]],df_tips = df[[tip_name_var]],tree=tr)

  # Order dataframe properly, if not the predictions will be wrong
  df<- df %>% .[match(tr$tip.label,.[[tip_name_var]]),]

  # Check if you want to use the best model
  if(model=="MF"){
    print("Performing model finding to use best fitting model")
    model = find_best_asr_model(df=df, tr=tr, tip_name_var=tip_name_var, pheno=pheno, node_states = node_states) %>% .$best_model
    print(paste0("Best model is: ",model))
  }

  # Run corHMM to estimate hidden rates
  corHMM_out = corHMM::corHMM(phy=tr, data=df[, c(tip_name_var,pheno)], rate.cat = 1, model=model, node.states = node_states,upper.bound=upper_bound,lower.bound=lower_bound)

  # Check if rates are at min or max bound
  check_rates_at_local_max(corHMM_out,upper_bound=upper_bound,lower_bound=lower_bound)

  # Model stats
  corHMM_model_summary = characterize_asr_model(corHMM_out)

  # Get Parent Child data
  outcome_str <- df[,pheno] %>% `names<-`(df[[tip_name_var]])
  parent_child_df <- get_parent_child_data(tr=tr, anc_data=corHMM_out$states, pheno_data=outcome_str, conf_threshold = conf_threshold, node_states=node_states)

  # Annotate parent child data
  parent_child_df <- get_phenotypic_continuation_data(parent_child_df,node_states=node_states)

  asr_output = list(corHMM_out = corHMM_out,corHMM_model_summary=corHMM_model_summary,parent_child_df = parent_child_df,node_states=node_states)
  return(asr_output)
}

check_phenotype <- function(phenotype_var){
  if(sum(unique(as.numeric(phenotype_var)) %in% c(0,1))!=2){
    stop("Phenotype is not formatted as binary variable with event = 1 and no-event = 0")
  }
}

check_tree <- function(tree){
  if(is.rooted(tree) == F){
    stop("Tree is not rooted. Please root your tree using functions like ape::root() or phytools::midpoint.root()")
  }
}

check_phenotype_tree_length <- function(pheno_var,df_tips,tree){
  if(length(pheno_var) != length(tree$tip.label)){
    stop("Tree length does not equal length of phenotype variable. Ensure the number of tips in the tree is the number of entries in your phenotype variable.")
  }
  if(sum(! tree$tip.label %in% df_tips)>0){
    stop("At least one tree tip name is not found in the dataset")
  }
}

check_rates_at_local_max <- function(corHMM_out,upper_bound,lower_bound){
  if(sum(corHMM_out$solution%in% as.character(c(upper_bound,lower_bound))) >0){
    warning("Rates are at the upper or lower maximum, consider updating the maximums to get a more accurate solution")
  }
}

get_parent_child_data <- function(tr, anc_data, pheno_data, conf_threshold=0.875, node_states){
  # Tree edge info
  edge <- tr$edge %>% as.data.frame
  colnames(edge) <- c("parent", "child")

  # Prediction data
  if(node_states == "marginal"){
    anc_data = as.data.frame(anc_data) %>% `colnames<-`(levels(as.factor(pheno_data)))
    anc_data[,"pred"] <- ifelse(anc_data[, 2] > conf_threshold, 1, ifelse(anc_data[, 1] >conf_threshold,0,0.5))
  } else if (node_states == "joint"){
    state_coding = as.factor(pheno_data) %>% levels %>% `names<-`(1:length(.))
    anc_data = data.frame(pred = anc_data)
    anc_data$pred =  recode(anc_data$pred, !!!as.list(state_coding)) %>% as.numeric
  } else {
    stop("Phylosuite is only optimized for marginal and joint corHMM")
  }

  # Add node data
  anc_data[,"node"] <-  1:nrow(anc_data)  + (nrow(anc_data)+1)
  num_edges <- nrow(edge)
  internal_nodes <- c(length(pheno_data) +1):(unlist(edge) %>% max)

  # Dataframe construction
  edge$parent_val <- rep(NA, num_edges)
  edge$child_val <- rep(NA, num_edges)

  # Prediction values for the parent & child
  for(i in internal_nodes){
    edge[edge$parent == i, "parent_val"] <- anc_data[anc_data$node == i, "pred"]
    edge[edge$child == i, "child_val"] <- anc_data[anc_data$node == i, "pred"]
  }

  # Internal nodes
  if(node_states == "marginal"){
    edge$child_mle <- rep(NA, num_edges)
    edge$parent_mle <- rep(NA ,num_edges)
    for(i in internal_nodes){
    edge[edge$parent == i, "parent_mle"] <- anc_data[anc_data$node == i, levels(as.factor(pheno_data))[2]] *100
    edge[edge$child == i, "child_mle"] <- anc_data[anc_data$node == i, levels(as.factor(pheno_data))[2]] *100
    }
  }

  # Values and name for the tip
  for(i in 1:(min(internal_nodes)-1)){
    edge[edge$child == i, "child_name"] <- tr$tip.label[[i]]
    edge[edge$child == i, "child_val"] <- subset(pheno_data,names(pheno_data)==edge[edge$child == i, "child_name"])
  }

  return(edge)
}

get_phenotypic_continuation_data <- function(parent_child_df,node_states){
  if(node_states =="joint"){
    parent_child_df <- parent_child_df %>% dplyr::mutate(transition = ifelse(parent_val != child_val, 1, 0),
                                                         gain = ifelse(child_val ==1 & parent_val==0, 1, 0),
                                                         loss = ifelse(child_val ==0 & parent_val==1, 1, 0),
                                                         continuation  =  ifelse(parent_val == child_val, 1, 0),
                                                         continuation_present = ifelse(continuation==1 & child_val==1,1,0),
                                                         continuation_absent = ifelse(continuation==1 & child_val==0,1,0))
  }
  if(node_states == "marginal"){
    parent_child_df <- parent_child_df %>% dplyr::mutate(transition = ifelse(parent_val != child_val, 1, 0),
                                                         transition_high = ifelse(parent_val ==0 & child_val==1 | parent_val==1 & child_val ==0, 1, 0),
                                                         transition_low = ifelse(parent_val ==0.5 & child_val==1 | parent_val==0.5 & child_val ==0, 1, 0),
                                                         gain = ifelse(child_val ==1 & parent_val==0.5 | child_val ==1 & parent_val==0, 1, 0),
                                                         gain_high  =  ifelse(child_val ==1 & parent_val==0, 1, 0),
                                                         gain_low = ifelse(child_val ==1 & parent_val==0.5, 1, 0),
                                                         loss = ifelse(child_val ==0 & parent_val==1 | child_val ==0 & parent_val==0.5, 1, 0),
                                                         loss_high  =  ifelse(child_val ==0 & parent_val==1, 1, 0),
                                                         loss_low  =  ifelse(child_val ==0 & parent_val==0.5, 1, 0),
                                                         continuation  =  ifelse(parent_val == child_val, 1, 0),
                                                         continuation_high  =  ifelse(parent_val == 1 & child_val ==1 | parent_val == 0 & child_val ==0,1,0),
                                                         continuation_low =  ifelse(parent_val == 0.5 & child_val ==0.5, 1, 0),
                                                         continuation_present = ifelse(continuation==1 & child_val==1,1,0),
                                                         continuation_present_high = ifelse(continuation_high==1 & child_val==1,1,0),
                                                         continuation_absent = ifelse(continuation==1 & child_val==0,1,0),
                                                         continuation_absent_high = ifelse(continuation_high==1 & child_val==0,1,0))
  }
  return(parent_child_df)
}
