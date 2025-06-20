# nearest_neighbor_algorithm: Nearest neighbor algorithm
#'
#' Identify the phylogenetic nearest neighbor of a tip
#'
#' This function uses phylogenetics to identify an isolate's nearest neighbor.
#' @param isolate Isolate of interest
#' @param tip_name_variable Name of variable containing tip names in the dataframe (metadata). Tip name variable must correspond to tip names in the tree.
#' @param phylogenetic_distance Matrix of phylogenetic distance
#' @param variant_distance Matrix with secondary distance measure, presumably SNPs. This is OPTIONAL. If not specified, the first entry will be gathered (OPTIONAL)
#' @param metadata Dataframe with trait and variables of interest
#' @param annotate Boolean (TRUE/FALSE) indicating if annotation is necessary
#' @param variables_of_interest Variables found in the metadata
#' @param comparison Boolean (TRUE/FALSE) indicating whether to focus on isolates with/without a trait
#' @param comparison_feature Name of variable to compare.
#' @return The nearest neighbor of interest.
#' @importFrom stringr str_split
#' @importFrom dplyr select
#' @export
nearest_neighbor_algorithm <- function(isolate, tip_name_variable, phylogenetic_distance, variant_distance, metadata, annotate = NULL, variables_of_interest = NULL, comparison = NULL, comparison_feature = NULL) {
  # Check entry values
  ## Check if nearest neighbor distance matrices are
  check_if_nearest_neighbor_data_is_matrix(phylogenetic_distance = phylogenetic_distance, variant_distance = variant_distance)
  ## Check if variables are in metadata
  if (annotate == TRUE) {
    check_if_variables_are_present(metadata = metadata, variables_of_interest = variables_of_interest, comparison_feature = comparison_feature)
  }

  # Get nearest neighbor
  nearest_neighbor <- get_nearest_neighbor(isolate = isolate, tip_name_variable = tip_name_variable, phylogenetic_distance = phylogenetic_distance, variant_distance = variant_distance, metadata = metadata, comparison = comparison, comparison_feature = comparison_feature)
  # Annotate Data
  if (annotate == TRUE) {
    nearest_neighbor <- annotate_nearest_neighbor(nearest_neighbor = nearest_neighbor, metadata = metadata, variables_of_interest = variables_of_interest)
  }
  return(nearest_neighbor)
}

# Get the nearest neighbor
get_nearest_neighbor <- function(isolate, tip_name_variable, phylogenetic_distance, variant_distance, metadata, comparison = NULL, comparison_feature = NULL) {
  # Get nearest neighbor candidates
  dist_mat <- phylogenetic_distance[rownames(phylogenetic_distance) != isolate, isolate]
  if (comparison == TRUE) {
    isolate_feature_value <- metadata[metadata[[tip_name_variable]] == isolate, comparison_feature]
    nearest_neighbor_candidates <- metadata[metadata[[comparison_feature]] != isolate_feature_value, tip_name_variable]
    dist_mat <- dist_mat[names(dist_mat) %in% nearest_neighbor_candidates]
  }
  nearest_neighbor <- names(which.min(dist_mat))
  multiple_neighbors <- ifelse(length(nearest_neighbor) > 1, TRUE, FALSE)
  # If multiple nearest neighbors, get closest using a secondary matrix
  if (multiple_neighbors == TRUE && is.null(variant_distance) == FALSE) {
    variant_dist <-  variant_distance[rownames(variant_distance) != isolate & rownames(variant_distance) %in% nearest_neighbor, ]
    nearest_neighbor <- names(which.min(variant_dist))
  }
  if (length(nearest_neighbor) > 1) {
    nearest_neighbor <- nearest_neighbor[[1]]
  }

  # Final data on earest neighbor
  ## Phylogenetic distance
  nearest_neighbor_phylogenetic_distance <- phylogenetic_distance[nearest_neighbor, isolate]
  data_row <- data.frame(tip_name = isolate, nearest_neighbor, nearest_neighbor_phylogenetic_distance)

  ## Variant distance
  if (is.null(variant_distance) == FALSE) {
    data_row$nearest_neighbor_variant_distance <- variant_distance[nearest_neighbor, isolate]
  }

  if (comparison == TRUE) {
    data_row[["comparison_feature"]] <- comparison_feature
    data_row[["isolate_comparison_value"]] <- isolate_feature_value
    data_row[["nearest_neighbor_comparison_value"]] <-  metadata[metadata[[tip_name_variable]] == nearest_neighbor, comparison_feature]
  }

  return(data_row)
}

# Annotate the nearest neighbor
annotate_nearest_neighbor <- function(nearest_neighbor, metadata, variables_of_interest) {
  nearest_neighbor_info <- nearest_neighbor[, !(names(nearest_neighbor) %in% c("tip_name", "nearest_neighbor"))]

  # Get metadata for isolate and its neighbor
  variables <- unique(variables_of_interest, comparison_feature)
  isolate_metadata <- metadata[metadata[[tip_name_variable]] == nearest_neighbor[["tip_name"]], ][variables]
  neighbor_metadata <-  metadata[metadata[[tip_name_variable]] == nearest_neighbor[["nearest_neighbor"]], ][variables]
  colnames(neighbor_metadata) <- paste0(colnames(neighbor_metadata), "_nn")

  # Construct metadata output
  nearest_neighbor <- data.frame(tip_name = nearest_neighbor[["tip_name"]],
                                 isolate_metadata,
                                 nearest_neighbor = nearest_neighbor[["nearest_neighbor"]],
                                 neighbor_metadata,
                                 nearest_neighbor_info)
  return(nearest_neighbor)
}
