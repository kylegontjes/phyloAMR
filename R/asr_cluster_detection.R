#' Identification of trait emergence and clustering of trait across the phylogeny using transition data
#'
#' This function identifies episodes of trait emergence and loss across a phylogenetic tree.
#'
#' @param df Dataframe with tip name variable and phenotype
#' @param tr Phylogenetic tree
#' @param tip_name_variable Name of variable containing tip names in df
#' @param patient_id Name of variable containing patient IDs, can be combined with faux_clusters option to factor into whether a cluster should have >1 patient. (Optional)
#' @param parent_child_df Parent child dataframe from asr() object
#' @param node_states Whether the reconstruction was "joint" or "marginal"
#' @param confidence Whether to use 'high' (i.e., 0 -> 1) or 'low' (i.e., any transition) confidence transitions when determining clustering with marginal ancestral state reconstruction results. If the confidence_threshold value in asr() was > 0.5, set confidence as 'low'. Otherwise, set confidence as 'high'.
#' @param simplify_faux_clusters Boolean (i.e., TRUE/FALSE), whether to collapse faux clusters (i.e., clusters where 1 patient contributes all isolates) as singletons without distinction (Optional)
#' @param simplify_revertant Boolean (i.e., TRUE/FALSE). Whether to collapse revertant episodes as isolates without the trait in the cleaned text string
#' @param collapse_cluster Boolean (i.e., TRUE/FALSE). Whether to create a variable that collapses cluster calls into one category
#' @return A tip-only dataframe with inferences on the ancestral history of these strains. Can be merged with parent_child_df from asr() if desired \describe{
#'     \item{asr_cluster}{Character string indicating cluster calls (cluster_[node]), singleton calls, traits without the feature (no feature), and revertant cases at the tip (revertant_tip) or clusters of revertants (revertant_cluster_[node]). If patient_id != NULL, additional calls may be provided where a cluster contains only one patient (cluster_[node]_1pt_only)}
#'     \item{patient_id}{Character string with the patient ID, if provided}
#'     \item{asr_cluster_renamed}{Character string where asr_cluster string was renamed as cluster [no. X], singletons, no feature. Clusters are ordered by presentation via ggtree. If simplify_revertant == TRUE, revertants are collapsed as 'No feature' }
#'     \item{asr_cluster_collapsed}{Character string where asr_cluster string was collapsed into cluster, singleton, no feature, and revertant. If simplify_revertant == TRUE, revertants are collapsed as 'No feature' }
#'   }
#' @importFrom dplyr case_when
#' @importFrom dplyr mutate
#' @importFrom dplyr first
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
  tip_parent_child_data <- subset(parent_child_df, child <= length(tr$tip.label))
  tip_parent_child_data[["tip_name"]] <- tip_parent_child_data[["child_name"]]
  tip_data_df <- suppressMessages(dplyr::left_join(tip_parent_child_data, clustering_data))
  rownames(tip_data_df) <- tip_data_df$tip_name

  return(tip_data_df)
}

# Get clustering data
get_clustering_data <- function(isolate, parent_child_df, tr, root_node, node_states, confidence) {
  tip_data <- parent_child_df[!is.na(parent_child_df$child_name) & parent_child_df$child_name  == isolate, ]

  # Classify tips with the trait
  if (tip_data[["child_value"]] == 1) {
    classification <- classify_tips_with_trait(node_data = tip_data, parent_child_df = parent_child_df, tr = tr, root_node = root_node, node_states = node_states, confidence = confidence)
    if (classification != "singleton") {
      classification <-  paste0("cluster_", classification)
    }
  }

  # Classify tips without the trait
  if (tip_data[["child_value"]] == 0) {
    classification <- classify_tips_without_trait(node_data = tip_data, parent_child_df = parent_child_df, tr = tr, root_node = root_node, node_states = node_states, confidence = confidence)
    if (!classification %in% c("no feature", "revertant_tip")) {
      classification <-  paste0("revertant_cluster_", classification)
    }
  }

  return(classification)
}

# Characterize tips with traits
classify_tips_with_trait <- function(node_data, parent_child_df, tr, root_node, node_states, confidence) {
  # Reclassify the singletons
  if ((node_states == "joint" && node_data[["gain"]] == 1) || (node_states == "marginal" && node_data[[paste0("gain_", confidence)]] == 1)) {
    previous <- get_parent_node(parental_node = node_data[["parent"]],parent_child_df =  parent_child_df)
    if ((node_states == "joint" && previous[["loss"]] == 1) || (node_states == "marginal" && previous[[paste0("loss_", confidence)]] == 1)) {
      classification <- get_transition_node(node = previous[['parent']], parent_child_df = parent_child_df, tr = tr, root_node = root_node, node_states = node_states, confidence = confidence)
    } else {
      classification <- "singleton"
    }
  } else {
    classification <- get_transition_node(node = node_data[['child']], parent_child_df = parent_child_df, tr = tr, root_node = root_node, node_states = node_states, confidence = confidence)
  }

  # Under marginal, consider episodes of ambiguity as singletons
  if (node_states == "marginal" && classification == "unsure") {
    classification <- "singleton"
  }

  return(classification)
}

# Characterize tips with traits
classify_tips_without_trait <- function(node_data, parent_child_df, tr, root_node, node_states, confidence) {
  # Reclassify the revertants at the tip
  if ((node_states == "joint" && node_data[["loss"]] == 1) || (node_states == "marginal" && node_data[[paste0("loss_", confidence)]] == 1)) {
    previous <- get_parent_node(parental_node =  node_data[["parent"]], parent_child_df = parent_child_df)
    if ((node_states == "joint" && previous[["gain"]] == 1) || (node_states == "marginal" && previous[[paste0("gain_", confidence)]] == 1)) {
      classification <- get_transition_node(node = previous[['parent']], parent_child_df = parent_child_df, tr = tr, root_node = root_node, node_states = node_states, confidence =  confidence)
    } else {
      classification <- "revertant_tip"
    }
  } else {
    classification <- get_transition_node(node = node_data[['child']], parent_child_df = parent_child_df, tr = tr, root_node = root_node, node_states = node_states, confidence =  confidence)
    # This is important, because our focus is on the emergence of a trait.
    # No reversion would be detected if there was a stretch of no trait all the way from tip to the root
    if (classification == "root") {
      classification <- "no feature"
    }
  }
  # Under marginal, consider episodes of ambiguity as having feature, but not revertant
  if (node_states == "marginal" && classification == "unsure") {
    classification <- "no feature"
  }

  return(classification)
}

## Utilitey function - get parent's node
get_parent_node <- function(parental_node, parent_child_df) {
  parent_child_df[parent_child_df$child == parental_node, ]
}

get_transition_node <- function(node, parent_child_df, tr, root_node, node_states, confidence) {
  # Get ancestral data
  ## Ancestors
  path_to_root <- ape::nodepath(tr, from = node, to = root_node)
  ancestors <- path_to_root[!path_to_root %in% c(node,root_node)]
  ## Ancestor parent child dataset
  ancestor_parent_child_df <-  parent_child_df[match(ancestors, parent_child_df$child), ]

  # Get node where a transition occurred
  if (node_states == "joint") {
    # Focusing on transitions, permitting analysis of gains and losses
    startvalue <-  parent_child_df[parent_child_df$child == node, "transition"]
    ancestor_values <-  ancestor_parent_child_df[["transition"]]
    # Classify individuals without transitions as trait being present at root (i.e., tip-to-root walk indicates no transitions)
    if (sum(ancestor_values) == 0) {
      node <- "root"
    } else {
      node <- ancestors[first(which(ancestor_values != startvalue))]
    }
  } else if (node_states == "marginal") {
    if (confidence == "high") {
      ancestor_values <-  ancestor_parent_child_df[["transition_high"]]
      startvalue <-   parent_child_df[parent_child_df$child == node, "transition_high"]
    } else {
      ancestor_values <-  ancestor_parent_child_df[["transition"]]
      startvalue <-  parent_child_df[parent_child_df$child == node, "transition"]
    }

    node <- ancestors[first(which(ancestor_values != startvalue))]

    if (is.na(node)) {
      node <- "unsure"
    }
  }
  return(node)
}

# Account for faux clusters
account_for_faux_clusters <- function(cluster_data) {
  # Determine cluster size
  cluster_calls <- subset(cluster_data$asr_cluster, grepl("cluster_", cluster_data$asr_cluster))
  cluster_size <- table(cluster_calls)
  # Isolates where gain at internal node was called, but no other isolate was connected to it
  true_clusters <- names(subset(cluster_size, cluster_size > 1))
  not_cluster <- names(cluster_size[!names(cluster_size) %in% true_clusters])
  # Clusters with only one patient
  cluster_size_pts <- sapply(true_clusters, FUN = function(x) {
    length(unique(cluster_data[cluster_data$asr_cluster == x, "patient_id"]))
  })
  not_more_than_one_pt <- names(subset(cluster_size_pts, cluster_size_pts == 1))
  # Relabel these as either (1) singletons/no feature or (2) clusters of 1 patient only
  ## Not clusters
  cleaned_cluster_string <- ifelse(cluster_data$asr_cluster %in% not_cluster, ifelse(grepl("revertant",cluster_data$asr_cluster), "no feature", "singleton"), cluster_data$asr_cluster)
  ## 1pt only clusters
  cleaned_cluster_string <- ifelse(cleaned_cluster_string %in% not_more_than_one_pt, paste0(cleaned_cluster_string, "_1pt_only"), cleaned_cluster_string)

  return(cleaned_cluster_string)
}

# Simplify clustering string
simplify_clustering_string <-  function(clustering_data, tr, simplify_faux_clusters, simplify_revertant) {
  # Get initial order
  initial_order <- clustering_data$tip_name
  #Reorder vector to tree plotting
  clustering_data <- clustering_data[match(ggtree::get_taxa_name(ggtree::ggtree(tr)), clustering_data$tip_name), c("tip_name", "asr_cluster")]

  # Convert Names to Easy to Report String
  clustering_data <- clustering_data %>% mutate(asr_cluster_renamed = case_when(asr_cluster == "no feature" ~ "No feature",
                                                                                asr_cluster == "singleton" ~ "Singleton",
                                                                                TRUE ~ asr_cluster))
  # Name clusters by position on tree
  cluster_isolates <- unique(subset(clustering_data$asr_cluster, grepl("cluster", clustering_data$asr_cluster) & !grepl("1pt|revertant", clustering_data$asr_cluster)))
  num_clusters <- length(cluster_isolates)

  # Name clusters
  if (num_clusters > 0) {
    clusters_renaming_vector <- create_renaming_vector(cluster_isolates, "Cluster")
    clustering_data$asr_cluster_renamed <- dplyr::recode(clustering_data$asr_cluster_renamed, !!! clusters_renaming_vector)
  }

  # Deal with faux clusters
  if (sum(grepl("1pt", clustering_data$asr_cluster)) > 0) {
    if (simplify_faux_clusters == FALSE) {
      faux_cluster_isolates <- unique(subset(clustering_data$asr_cluster, grepl("1pt", clustering_data$asr_cluster)))
      faux_cluster_renaming_vector <- create_renaming_vector(faux_cluster_isolates, "Single patient cluster")
      clustering_data$asr_cluster_renamed <- dplyr::recode(clustering_data$asr_cluster_renamed, !!! faux_cluster_renaming_vector)
    } else {
      clustering_data$asr_cluster_renamed <- ifelse(grepl("1pt", clustering_data$asr_cluster_renamed), "Singleton", clustering_data$asr_cluster_renamed)
  }
  }

  # Deal with revertants
  if (sum(grepl("revertant", clustering_data$asr_cluster)) > 0) {
    if (simplify_revertant == TRUE) {
      # Convert all revertant episodes to 'no feature'
      clustering_data <- clustering_data %>% mutate(asr_cluster_renamed = case_when(grepl("revertant_", asr_cluster) ~ "No feature", TRUE ~ asr_cluster_renamed))
    } else {
      # Rename revertant_tip
      clustering_data <- clustering_data %>% mutate(asr_cluster_renamed = case_when(asr_cluster == "revertant_tip" ~ "Revertant tip", TRUE ~ asr_cluster_renamed))
    # Rename revertant clusters
    if (sum(grepl("revertant_", clustering_data$asr_cluster) & !grepl("revertant_tip", clustering_data$asr_cluster)) > 0) {
      revertant_cluster_isolates <-  unique(subset(clustering_data$asr_cluster, grepl("revertant_", clustering_data$asr_cluster) & !grepl("revertant_tip", clustering_data$asr_cluster)))
      revertant_cluster_renaming_vector <- create_renaming_vector(revertant_cluster_isolates, "Revertant Cluster")
      clustering_data$asr_cluster_renamed <- dplyr::recode(clustering_data$asr_cluster_renamed, !!! revertant_cluster_renaming_vector)
    }
  }
  }

  # Match order back to original dataframe
  clustering_data <- clustering_data[match(initial_order, clustering_data$tip_name), ]
  return(clustering_data$asr_cluster_renamed)
}

## Utility function: create renaming vector for simplify_clustering_string
create_renaming_vector <- function(string, naming_scheme) {
  renamed_string <- paste0(naming_scheme, " ", seq_len(length(string)))
  names(renamed_string) <- string
  return(renamed_string)
}

# Group by category
group_by_category <- function(string, simplify_faux_clusters, simplify_revertant) {
  collapsed_string <- case_when(grepl("Cluster", string) & !grepl("Revertant|Single patient", string) ~ "Cluster",
                                TRUE ~ string)

  # Simplify revertants will collapse values as 'No feature'
  if (simplify_revertant == FALSE) {
    collapsed_string <- case_when(grepl("Revertant", collapsed_string) ~ "Revertant",
                                  TRUE ~ collapsed_string)
  } else {
    collapsed_string <- case_when(grepl("Revertant", collapsed_string) ~ "No feature",
                                  TRUE ~ collapsed_string)
  }

  # If simplify_faux_clusters was false, these will still need to be grouped as singletons
  if (simplify_faux_clusters == FALSE) {
    collapsed_string <- case_when(grepl("Single patient",  collapsed_string) ~ "Singleton",
                                  TRUE ~ collapsed_string)
  }
  return(collapsed_string)
}

