#' Cluster detection
#'
#' Description of what the function does.
#'
#' @param df Dataframe with tip name variable and phenotype
#' @param tr Phylogenetic tree
#' @param tip_name_variable Name of variable containing tip names in df
#' @param patient_id Name of variable containing patient IDs, can be combined with faux_clusters option to factor into whether a cluster should have >1 patient. (Optional)
#' @param parent_child_df Parent child dataframe from asr() object
#' @param node_states Whether the reconstruction was "joint" or "marginal"
#' @param confidence Whether to use 'high' (i.e., 0 -> 1) or 'low' (0 -> 0.5) confidence transitions when determining clustering. ONLY USAGE FOR MARGINAL STATE RECONSTRUCTIONS
#' @param simplify_faux_clusters Booleane (i.e., TRUE/FALSE), whether to collapse faux clusters (i.e., clusters where 1 patient contributes all isolates) as singletons without distinction (Optional)
#' @param simplify_revertant Boolean (i.e., TRUE/FALSE). Whether to collapse revertant episodes as isolates without the trait in the cleaned text string
#' @param collapse_cluster Boolean (i.e., TRUE/FALSE). Whether to create a variable that collapses cluster calls into one category
#' @return A tip-only dataframe with inferences on the history of these strains. Can be merged with parent_child_df from asr() if desired
#' @export
asr_cluster_detection <- function(df, tr, tip_name_variable, patient_id = NULL, parent_child_df, node_states = "joint", confidence = NULL, simplify_faux_clusters = FALSE, simplify_revertant = TRUE, collapse_cluster = TRUE) {
  # Check if states are as desired
  check_joint_confidence(node_states, confidence)

  # Check faux_cluster
  check_faux_clusters(patient_id, simplify_faux_clusters)

  # Set up data frame
  clustering_data <- data.frame(tip_name = df[[tip_name_variable]])
  rownames(clustering_data) <- clustering_data$tip_name

  # Get clustering data
  root_node <- length(tr$tip.label) + 1
  clustering_data$asr_cluster <- sapply(clustering_data$tip_name, FUN = get_clustering_data, parent_child_df = parent_child_df, tr = tr, root_node = root_node, node_states = node_states, confidence = confidence)

  # If patient ID is not null, account for faux clusters (e.g., relabel clusters with only one patient)
  if (!is.null(patient_id)) {
    clustering_data$patient_id <- df[[patient_id]]
    clustering_data$asr_cluster <- account_for_faux_clusters(clustering_data)
  }

  # Simplify clustering output
  clustering_data$asr_cluster_renamed <- simplify_clustering_string(clustering_data = clustering_data, tr = tr, simplify_faux_clusters = simplify_faux_clusters, simplify_revertant = simplify_revertant)

  # Collapse cluster output through grouping by category
  if (collapse_cluster == TRUE) {
    clustering_data$asr_cluster_collapsed <- group_by_category(string = clustering_data$asr_cluster_renamed, simplify_faux_clusters = simplify_faux_clusters, simplify_revertant = simplify_revertant)
  }

  # Get tip data and merge it with the clustering data
  tip_data <- subset(parent_child_df, child <= ape::Nnode(tr) + 1)
  tip_data[["tip_name"]] <- tip_data[["child_name"]]
  tip_data_df <- suppressMessages(dplyr::left_join(tip_data, clustering_data))
  rownames(tip_data_df) <- tip_data_df$tip_name

  return(tip_data_df)
}

# Get clustering data
get_clustering_data <- function(isolate, parent_child_df, tr, root_node, node_states, confidence) {
  tip_row <- parent_child_df[!is.na(parent_child_df$child_name) & parent_child_df$child_name  == isolate, ]

  # Classify tips with the trait
  if (tip_row[['child_val']] == 1) {
    classification <- classify_tips_with_trait(node_data = tip_row, parent_child_df = parent_child_df, tr = tr, root_node = root_node, node_states = node_states, confidence = confidence)
    if(classification != 'singleton'){
      classification <-  paste0("cluster_",classification)
    }
  }

  # Classify tips without the trait
  if(tip_row[['child_val']] == 0){
    classification <- classify_tips_without_trait(node_data = tip_row, parent_child_df = parent_child_df, tr = tr, root_node = root_node, node_states = node_states, confidence = confidence)
    if (!classification %in% c("no feature", "revertant_tip")) {
      classification <-  paste0("revertant_cluster_",classification)
    }
  }

  return(classification)
}

# Characterize tips with traits
classify_tips_with_trait <- function(node_data, parent_child_df, tr, root_node, node_states, confidence) {
  # Reclassify the singletons
  if ((node_states =="joint" && node_data[['gain']] == 1) || (node_states =="marginal" && node_data[[paste0('gain_',confidence)]] == 1)){
    previous <- get_parent_node(node_data[['parent']], parent_child_df)
    if ((node_states == "joint" && previous[['loss']] ==1) || (node_states == "marginal" && previous[[paste0('loss_',confidence)]] == 1)){
      classification <- get_transition_node(previous, parent_child_df,tr, root_node = root_node, node_states, confidence)
    } else {
      classification <- "singleton"
    }
  } else{
    classification <- get_transition_node( node_data, parent_child_df,tr, root_node = root_node, node_states, confidence)
  }
  return(classification)
}

# Characterize tips with traits
classify_tips_without_trait <- function(node_data, parent_child_df, tr, root_node, node_states, confidence) {
  # Reclassify the revertants at the tip
  if ((node_states =="joint" && node_data[['loss']] == 1) || (node_states =="marginal" && node_data[[paste0('loss_',confidence)]] == 1)){
    previous <- get_parent_node(node_data[['parent']], parent_child_df)
    if ((node_states == "joint" && previous[['gain']] ==1 || node_states == "marginal" && previous[[paste0('gain_',confidence)]] == 1)){
      classification <- get_transition_node(node_data, parent_child_df, tr, root_node = root_node, confidence)
    } else {
      classification <- "revertant_tip"
    }
  } else{
    classification <- get_transition_node(node_data, parent_child_df, tr,root_node = root_node, node_states, confidence)
    # This is important, because our focus is on the emergence of a trait.
    # No reversion would be detected if there was a stretch of no trait all the way from tip to the root
    if(classification == 'root'){
      classification <- 'no feature'
    }
  }
  return(classification)
}

## Utilitey function - get parent's node
get_parent_node <- function(parental_node, parent_child_df) {
  parent_child_df[parent_child_df$child == parental_node,]
}

# Account for faux clusters
account_for_faux_clusters <- function(cluster_data) {
  # Determine cluster size
  cluster_calls <- subset(cluster_data$asr_cluster, grepl("cluster_", cluster_data$asr_cluster))
  cluster_size <- table(cluster_calls)
  # Isolates where gain at internal node was called, but no other isolate was connected to it
  true_clusters <- names(subset(cluster_size, cluster_size>1))
  not_cluster <- names(cluster_size[!names(cluster_size) %in% true_clusters])
  # Clusters with only one patient
  cluster_size_pts <- sapply(true_clusters, FUN = function(x) {
    length(unique(cluster_data[cluster_data$asr_cluster == x,'patient_id']))
  })
  not_more_than_one_pt <- names(subset(cluster_size_pts, cluster_size_pts == 1))
  # Relabel these as either (1) singletons or (2) clusters of 1 patient only
  cleaned_cluster_string <- case_when(
    cluster_data$asr_cluster %in% not_cluster ~ "singleton",
    cluster_data$asr_cluster %in% not_more_than_one_pt ~ paste0(cluster_data$asr_cluster, "_1pt_only"),
    TRUE ~ cluster_data$asr_cluster)
  return(cleaned_cluster_string)
}

# Simplify clustering string
simplify_clustering_string <-  function(clustering_data, tr, simplify_faux_clusters, simplify_revertant) {
  # Get initial order
  initial_order <- clustering_data$tip_name
  #Reorder vector to tree plotting
  clustering_data <- clustering_data[match(ggtree::get_taxa_name(ggtree::ggtree(tr)), clustering_data$tip_name), c('tip_name','asr_cluster')]

  # Convert Names to Easy to Report String
  clustering_data <- clustering_data %>% mutate(asr_cluster_renamed = case_when(asr_cluster == 'no feature' ~ "No feature",
                                                                                asr_cluster == 'singleton' ~ 'Singleton',
                                                                                TRUE ~ asr_cluster))
  # Name clusters by position on tree
  cluster_isolates <- unique(subset(clustering_data$asr_cluster,grepl("cluster",clustering_data$asr_cluster) & !grepl("1pt|revertant",clustering_data$asr_cluster)))
  num_clusters <- length(cluster_isolates)

  # Name clusters
  if (num_clusters >0){
    clusters_renaming_vector <- create_renaming_vector(cluster_isolates,'Cluster')
    clustering_data$asr_cluster_renamed <- dplyr::recode(clustering_data$asr_cluster_renamed, !!! clusters_renaming_vector)
  }

  # Deal with faux clusters
  if (sum(grepl("1pt",clustering_data$asr_cluster)) > 0) {
  if (simplify_faux_clusters == FALSE) {
    faux_cluster_isolates <- unique(subset(clustering_data$asr_cluster,grepl("1pt",clustering_data$asr_cluster)))
    faux_cluster_renaming_vector <- create_renaming_vector(faux_cluster_isolates,'Single patient cluster')
    clustering_data$asr_cluster_renamed <- dplyr::recode(clustering_data$asr_cluster_renamed, !!! faux_cluster_renaming_vector)
  } else {
    clustering_data$asr_cluster_renamed <- ifelse(grepl("1pt", clustering_data$asr_cluster_renamed), "Singleton", clustering_data$asr_cluster_renamed)
  }
  }

  # Deal with revertants
  if (sum(grepl("revertant",clustering_data$asr_cluster)) > 0) {
  if (simplify_revertant == TRUE) {
    # Convert all revertant episodes to 'no feature'
    clustering_data <- clustering_data %>% mutate(asr_cluster_renamed = case_when(grepl("revertant_",asr_cluster) ~ "No feature", TRUE ~ asr_cluster_renamed))
  } else {
    # Rename revertant_tip
    clustering_data <- clustering_data %>% mutate(asr_cluster_renamed = case_when(asr_cluster == "revertant_tip" ~ "Revertant tip", TRUE ~ asr_cluster_renamed))
    # Rename revertant clusters
    if (sum(grepl("revertant_",clustering_data$asr_cluster) & !grepl("revertant_tip",clustering_data$asr_cluster)) > 0) {
      revertant_cluster_isolates <-  unique(subset(clustering_data$asr_cluster,grepl("revertant_",clustering_data$asr_cluster) & !grepl("revertant_tip",clustering_data$asr_cluster)))
      revertant_cluster_renaming_vector <- create_renaming_vector(revertant_cluster_isolates,'Revertant Cluster')
      clustering_data$asr_cluster_renamed <- dplyr::recode(clustering_data$asr_cluster_renamed, !!! revertant_cluster_renaming_vector)
    }
  }
  }

  # Match order back to original dataframe
  clustering_data <- clustering_data[match(initial_order, clustering_data$tip_name), ]
  return(clustering_data$asr_cluster_renamed)
}

## Utility function: create renaming vector for simplify_clustering_string
create_renaming_vector <- function(string,naming_scheme){
  renamed_string <- paste0(naming_scheme, " ", seq_len(length(string)))
  names(renamed_string) <- string
  return(renamed_string)
}

# Group by category
group_by_category <- function(string, simplify_faux_clusters, simplify_revertant) {
  collapsed_string <- case_when(grepl("Cluster",string) & !grepl("Revertant|Single patient", string) ~ "Cluster",
                                TRUE ~ string)

  if (simplify_revertant == FALSE){
    collapsed_string <- case_when(grepl("Revertant", collapsed_string) ~ "Revertant",
                                  TRUE ~ collapsed_string)
  } else {
    collapsed_string <- case_when(grepl("Revertant", collapsed_string) ~ "No feature",
                                  TRUE ~ collapsed_string)
  }

  if(simplify_faux_clusters == FALSE){
    collapsed_string <- case_when(grepl("Single patient",  collapsed_string) ~ "Singleton",
                                  TRUE ~ collapsed_string)
  }
  return(collapsed_string)
}
get_transition_node <- function(child_data, parent_child_df, tr, root_node, node_states, confidence){
  # Get ancestral data
  ancestors <- nodepath(tr,from = child_data$child, to = root_node)
  ancestors <- ancestors[ancestors != root_node]
  ancestor_vals <-  parent_child_df[match(ancestors, parent_child_df$child,),'child_val']

  # Get stretches to the root
  startvalue = child_data$child_val
  if(node_states == 'joint'){
    end = ancestors[last(which(ancestor_vals == startvalue))]
    if(end == last(ancestors)){
      end = "root"
    }
  }
  return(end)
}
