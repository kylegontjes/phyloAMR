### Checks ###
#### Ancestral state reconstruction checks ####
## Check phenotype is formatted
check_trait <- function(trait_var) {
  if (sum(unique(as.numeric(trait_var)) %in% c(0, 1)) != 2) {
    stop("Phenotype is not formatted as binary variable with event = 1 and no-event = 0")
  }
}

## Check tree has relevant contents
check_tree <- function(tree) {
  if (ape::is.rooted(tree) == FALSE) {
    stop("Tree is not rooted. Please root your tree using functions like ape::root() or phytools::midpoint.root()")
  }
}

## Check length of phenotype matches number of tips in tree
check_trait_tree_length <- function(trait_var, df_tips, tree) {
  if (length(trait_var) != length(tree$tip.label)) {
    stop("Tree length does not equal length of phenotype/trait variable. Ensure the number of tips in the tree is the number of entries in your phenotype/trait variable.")
  }
  if (sum(! tree$tip.label %in% df_tips) > 0) {
    stop("At least one tree tip name is not found in the dataset")
  }
}

## Joint and confidence
### Check joint and confidence value for asr transition statistics
check_joint_confidence_value <- function(node_states, confidence_threshold){
  if(node_states == "joint" & is.null(confidence_threshold)== FALSE) {
    stop("When using joint reconstruction states, you must specify confidence as NULL")
  }
  if(node_states == 'marginal'){
    if(is.numeric(confidence_threshold) == FALSE){
      stop("Confidence threshold must be numeric and between 1 and 0")
    }
    if(confidence_threshold >1 | confidence_threshold < 0) {
      stop("Confidence threshold must be between 1 and 0")
  }
  }
}

### Check joint confidence for asr_cluster_detection calling
check_joint_confidence <- function(node_states, confidence){
  if(node_states == "joint" & is.null(confidence) == FALSE) {
    stop("When using joint reconstruction states, you must specify confidence as NULL")
  }
}

## Check if the corHMM rates reached upper or lower bounds
check_rates_at_local_max <- function(corHMM_output, upper_bound, lower_bound) {
  if (sum(corHMM_output$solution %in% as.character(c(upper_bound, lower_bound))) > 0) {
    warning("Rates are at the upper or lower maximum, consider updating the maximums to get a more accurate solution")
  }
}

#### ASR transitional data analysis ####
## Check if appropriate columns exist when requesting statistics for marginal reconstruction
check_marginal_columns <- function(node_states, parent_child_df) {
  if (node_states == "marginal" && c(!"transition_high" %in% colnames(parent_child_df))) {
    stop("Does not have marginal columns, ensure your parent_child_df is constructed using marginal reconstruction")
  }
  if (node_states == "joint" && "transition_high" %in% colnames(parent_child_df)) {
    stop("This is a marginal reconstruction, specify node_state as 'marginal'")
  }
}

#### Cluster detection analyses ####
## Faux clusters requiring patient_id information
check_faux_clusters <- function(patient_id, simplify_faux_clusters){
  if(is.null(patient_id) == TRUE & simplify_faux_clusters != FALSE) {
    stop("Must have patient id variable for when faux clusters is specified")
  }
}
