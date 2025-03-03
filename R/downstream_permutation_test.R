#' downstream_permutation_test: Permutation test for downstream analysis
#'
#' Permutation test for analysis of downstream gain and loss of a comparitor trait on stretches of a different trait of interest
#'
#' @param comparitor Comparitor trait, such as a genotype
#' @param df Dataframe with trait and comparitor and tip_name_var
#' @param tr Phylogenetic tree
#' @param tip_name_var Tip name variable
#' @param node_states Joint or marginal reconstruction
#' @param conf_threshold Set value for confidence threshold of MLE calls if marginal
#' @param confidence Whether to use high or low confidence transition nodes when node_states are marginal
#' @return Summary stats for downstream gain and loss of a trait
#' @export
downstream_permutation_test <- function(comparitor,df,tr,tip_name_var,trait,node_states = "joint",conf_threshold=NULL,confidence=NULL,num_permutations=1000,num_cores = 6){
  library(pbmcapply)
  # trait
  trait_asr = asr(df=df,tr=tr,tip_name_var=tip_name_var,pheno = trait,model = "ER",node_states = node_states,conf_threshold = conf_threshold) %>% .$parent_child_df

  # comparitor
  comparitor_asr = asr(df=df,tr=tr,tip_name_var=tip_name_var,pheno = comparitor,model = "ER",node_states = node_states,conf_threshold = conf_threshold) %>% .$parent_child_df

  # Get observed
  trait_downstream_obs = downstream_gain_loss(comparitor_asr,trait_asr,tr,node_states = node_states,confidence=confidence)
  # permute
  num_isolates = nrow(df)
  num_features = sum(as.numeric(df[[comparitor]]))

  trait_runs <- pbmclapply(seq_len(num_permutations),FUN=function(x){sample(df[[comparitor]],num_isolates,replace=F)},mc.cores=num_cores) %>% do.call(cbind,.) %>% data.frame
  rownames(trait_runs) = df[,tip_name_var]
  permutation_names = paste0("p",seq_len(num_permutations))
  colnames(trait_runs) = permutation_names
  trait_runs=data.frame(trait_runs,tip_name_var = df[,tip_name_var])
  asr_permutation = pbmclapply(permutation_names,FUN=function(x){
    asr_result = asr(df=trait_runs,tr=tr,tip_name_var="tip_name_var",pheno = x,model = "ER",node_states = node_states,conf_threshold = confidence) %>% .$parent_child_df
    downstream_perm = downstream_gain_loss(asr_result,trait_asr,tr)

    return(downstream_perm)
  },mc.cores = num_cores) %>% do.call(rbind,.)

  trait_downstream_obs$gain_pval = c(1+sum(asr_permutation$gains_num >= trait_downstream_obs$gains_num)) / c(1+num_permutations)
  trait_downstream_obs$gains_strecthes_pval = c(1+sum(asr_permutation$stretches_w_gains_num >= trait_downstream_obs$stretches_w_gains_num)) / c(1+num_permutations)
  trait_downstream_obs$loss_pval = c(1+sum(asr_permutation$loss_num >= trait_downstream_obs$loss_num)) / c(1+num_permutations)
  trait_downstream_obs$loss_streches_pval = c(1+sum(asr_permutation$stretches_w_losses_num >= trait_downstream_obs$stretches_w_losses_num)) / c(1+num_permutations)

  return(trait_downstream_obs)
}
