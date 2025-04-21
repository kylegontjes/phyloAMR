# Nearest neighbor algorithm

#' Generate nearest neighbor using phylogenetic distance
#'
#' This function will use the phylogenetic tree to identify an isolate's nearest neighbor. This genetically
#' similar strain can be used later.

nearest_neighbor_algorithm <- function(isolate,phylogenetic_distance,variant_distance,metadata,annotate=NULL,variables_of_interest=NULL,comparison=NULL,comparison_feature=NULL){
  # Get nearest neighbor
  nearest_neighbor <- get_nearest_neighbor(isolate = isolate,phylogenetic_distance = phylogenetic_distance,variant_distance=variant_distance,metadata=metadata,comparison=comparison,comparison_feature=comparison_feature)
  # Annotate Data
  if(annotate==T){
    nearest_neighbor <- annotate_nearest_neighbor(isolate,nearest_neighbor=nearest_neighbor$nearest_neighbor,metadata=metadata,variables_of_interest=variables_of_interest,comparison_feature=comparison_feature)
  }
  return(nearest_neighbor)
}

get_nearest_neighbor <- function(isolate,phylogenetic_distance,variant_distance,metadata=metadata,comparison=NULL,comparison_feature=NULL){
  # Get nearest neighbor candidates
  dist_mat <- phylogenetic_distance[,isolate] %>% subset(names(.) != isolate)
  if(comparison == T){
    isolate_feature_value <- metadata %>% subset(isolate_no == isolate) %>% .[,comparison_feature]
    nn_candidates <- subset(metadata,get(comparison_feature) != isolate_feature_value) %>% .$isolate_no
    dist_mat <-  subset(dist_mat,names(dist_mat) %in% nn_candidates)
  }
  nearest_neighbor <- subset(dist_mat,dist_mat == min(dist_mat)) %>% names %>% paste0(collapse=",")
  multiple_neighbors <- ifelse(length(unlist(stringr::str_split(nearest_neighbor,",")))>1,T,F)
  # If multiple nearest neighbors, get closest using a secondary matrix
  if(multiple_neighbors == T){
    nn_candidates <- c(str_split(nearest_neighbor,",") %>% unlist)
    variant_mat <-  variant_distance %>% select(isolate) %>% subset(rownames(.) != isolate & rownames(.) %in% nn_candidates)
    nearest_neighbor <- subset(variant_mat,variant_mat == min(variant_mat)) %>% rownames %>% paste0(collapse=",")
  }
  if(length(unlist(stringr::str_split(nearest_neighbor,",")))>1){
    nearest_neighbor <- stringr::str_split(nearest_neighbor,",") %>% unlist %>% .[1]
  }
  data_row <- data.frame(isolate_no = isolate,nearest_neighbor=nearest_neighbor)
  return(data_row)
}

annotate_nearest_neighbor <- function(isolate,nearest_neighbor,metadata,variables_of_interest,comparison_feature){
  isolate_metadata <- metadata %>% subset(rownames(.) == isolate) %>% dplyr::select(dplyr::any_of(c("isolate_no",variables_of_interest,comparison_feature)))
  neighbor_metadata <- metadata %>% subset(rownames(.) == nearest_neighbor) %>% dplyr::select(variables_of_interest,comparison_feature) %>% `colnames<-`(paste0(colnames(.),"_nn"))
  nearest_neighbor <- data.frame(isolate_metadata,nearest_neighbor,neighbor_metadata)
  return(nearest_neighbor)
}
