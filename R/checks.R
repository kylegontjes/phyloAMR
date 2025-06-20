##### Checks #####
#### Ancestral state reconstruction ####
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

#### asr_cluster_detection ####
## Check joint confidence for asr_cluster_detection calling
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
    stop("Does not have marginal columns, ensure that your parent_child_df is constructed using marginal reconstruction")
  }
  if (node_states == "joint" && "transition_high" %in% colnames(parent_child_df)) {
    stop("This is a marginal reconstruction, specify node_states as 'marginal'")
  }
}

#### Cluster detection analyses ####
## Faux clusters requiring patient_id information
check_faux_clusters <- function(patient_id, simplify_faux_clusters){
  if(is.null(patient_id) == TRUE & simplify_faux_clusters != FALSE) {
    stop("Must have patient_id variable when simplify_faux_clusters is requested")
  }
}

#### asr cluster analysis ####
check_asr_content <- function(tip_data_df){
  if (sum(grepl("asr_cluster",colnames(tip_data_df))) == 0) {
    stop("The input does not have asr_cluster results. Double check input is the tip_data_df from asr_cluster_detection()")
  }
}

#### Phyloaware dataset curation check ####
### Check has asr and trait
check_asr_trait <- function(df,trait){
  if (sum(c(trait, "asr_cluster") %in% colnames(df)) != 2) {
    stop("The input does not have asr_cluster and/or trait. Check content.")
  }
}


check_patient_id_and_culture_date <- function(df, patient_id, culture_date){
  if(sum(c(patient_id, culture_date) %in% colnames(df)) != 2) {
    stop("The input does not have the provided patient_id and/or culture date variable. These are required when specifying first_present == TRUE")
  }
}

#### Nearest neighbor algorithm ####
check_if_nearest_neighbor_data_is_matrix <- function(phylogenetic_distance, variant_distance) {
  if (is.matrix(phylogenetic_distance) == FALSE) {
    stop("The phylogenetic distance object must be formatted as a matrix")
  }

  if (is.null(variant_distance) == FALSE & is.matrix(variant_distance) == FALSE) {
    stop("The variant distance object must be formatted as a matrix")
  }
}

check_if_variables_are_present <- function(metadata = metadata, variables_of_interest = variables_of_interest, comparison_feature = comparison_feature){
  if (sum(!c(variables_of_interest, comparison_feature) %in% colnames(metadata)) > 0) {
    stop("At least one requested variable was not included in the metadata object")
  }
}
