#' Analysis of downstream transitions of comparitor traits
#'
#' Downstream gain and loss of a comparitor trait on stretches of a different trait of interest
#'
#' @param comparitor_parent_child_df Parent child dataset for a comparitor trait, such as a genotype or a different trait/phenotype
#' @param trait_parent_child_df Parent child dataset for a trait of interest, such as a trait/phenotype
#' @param tr Phylogenetic tree
#' @param node_states Joint or marginal reconstruction
#' @param confidence Whether to use high or low confidence transition nodes when node_states are marginal
#' @return Summary stats for downstream gain and loss of a trait
#' @importFrom utils head
#' @export
downstream_transitions <- function(comparitor_parent_child_df, trait_parent_child_df, tr, node_states = "joint", confidence = NULL) {
  # Get downstream nodes from the gain eventss
  downstream_nodes <- get_trait_traces_on_tree(parent_child_df = trait_parent_child_df, tr = tr, node_states = node_states)
  # Get gain/losses on traces on the tree with the trait
  downstream_changes <- get_gain_loss_on_stretches(comparitor_parent_child_df = comparitor_parent_child_df, downstream_nodes = downstream_nodes, node_states = node_states, confidence = confidence)
  return(downstream_changes)
}

get_trait_traces_on_tree <- function(parent_child_df, tr, node_states) {
  all_possible_paths <- ape::nodepath(tr)

  # Get gain location
  if (node_states == "joint") {
    gains <- parent_child_df[parent_child_df$gain == 1, "child"]
  } else {
    confidence_variable <- ifelse(confidnece == 'high', "gain_high", "gain")
    gains <- parent_child_df[parent_child_df[[confidence_variable]] == 1, "child"]
  }

  # Decipher tip and ancestral gains
  ancestral_gains <- subset(gains, gains >= min(parent_child_df$parent))

  # Trait positions
  trait_found <- unique(c(parent_child_df[parent_child_df$child_val ==1, "child"], parent_child_df[parent_child_df$parent_val ==1, "parent"]))

  # Decipher node paths with gain events
  downstream_nodes <- sapply(ancestral_gains, FUN = function(gain, all_possible_paths){
    # Paths with the gain
    paths_w_gain <- all_possible_paths[sapply(all_possible_paths, FUN = function(path){gain %in% path})]

    # Paths from the gain
    gain_stretches <- lapply(paths_w_gain, FUN = function(gain_path){
      gain_path[which(gain_path == gain):length(gain_path)]
    })

    # Annotate with trait
    downstream_stretches_w_trait <- lapply(gain_stretches, FUN = function(stretch, trait_found){
      positive_nodes <- stretch %in% trait_found
      if(sum(positive_nodes) != length(stretch)){
        positive_stretch <- stretch[2:c(dplyr::first(which(positive_nodes != TRUE))-1)]
      } else{
        positive_stretch <- stretch[2:dplyr::last(which(positive_nodes == TRUE))]
      }
      return(positive_stretch)
    }, trait_found = trait_found)

    unique_gain_stretches <- unique(downstream_stretches_w_trait)

    # Get total downstream
    downstream_nodes <- unique(unlist(unique_gain_stretches))

    return(downstream_nodes)

  }, all_possible_paths = all_possible_paths)

  names(downstream_nodes) <- ancestral_gains

  return(downstream_nodes)
}

get_gain_loss_on_stretches <- function(comparitor_parent_child_df, downstream_nodes, node_states, confidence) {
  # Get gain location
  if (node_states == "joint") {
    comparitor_gains <- comparitor_parent_child_df[comparitor_parent_child_df$gain == 1, "child"]
    comparitor_losses <- comparitor_parent_child_df[comparitor_parent_child_df$loss == 1, "child"]
  } else {
    confidence_gain <- ifelse(confidnece == 'high', "gain_high", "gain")
    confidence_loss <- ifelse(confidnece == 'high', "loss_high", "loss")
    comparitor_gains <- comparitor_parent_child_df[comparitor_parent_child_df[[confidence_gain]] == 1, "child"]
    comparitor_losses <- comparitor_parent_child_df[comparitor_parent_child_df[[confidence_loss]] == 1, "child"]
  }

  # Get downstream positions
  downstream_positions <- unlist(downstream_nodes,use.names = FALSE)
  num_downstream_edges <- length(downstream_positions)

  # Downstream gain events
  gains <- sort(as.numeric(downstream_positions[which(downstream_positions %in% comparitor_gains)]))
  gains_str <- paste0(gains, collapse = ",")
  gains_num <- length(gains)
  gains_prop <- gains_num / num_downstream_edges

  # Downstream loss events
  loss <- sort(as.numeric(downstream_positions[which(downstream_positions %in% comparitor_losses)]))
  loss_str <- paste0(loss, collapse = ",")
  loss_num <- length(loss)
  loss_prop <- loss_num / num_downstream_edges

  # Number of paths
  num_stretches <- length(downstream_nodes)

  # Gains
  stretches_and_gains <- lapply(downstream_nodes, FUN = function(stretch, comparitor_gains){
    num_gains <- sum(stretch %in% comparitor_gains)
    ifelse(num_gains > 0, TRUE, FALSE)
  }, comparitor_gains = comparitor_gains)
  stretches_w_gains <- names(subset(stretches_and_gains, stretches_and_gains == TRUE))

  stretches_w_gains_str <- ifelse(length(stretches_w_gains) > 0, paste0(stretches_w_gains, collapse = ","), character(0))
  stretches_w_gains_num <- length(stretches_w_gains)
  stretches_w_gains_prop <- round(stretches_w_gains_num / num_stretches, 2)

  # Losses
  stretches_and_losses <- lapply(downstream_nodes, FUN = function(stretch, comparitor_losses){
    num_gains <- sum(stretch %in% comparitor_losses)
    ifelse(num_gains > 0, TRUE, FALSE)
  }, comparitor_losses = comparitor_losses)
  stretches_w_losses <- names(subset(stretches_and_losses, stretches_and_losses == TRUE))

  stretches_w_losses_str <- ifelse(length(stretches_w_losses) > 0, paste0(stretches_w_losses, collapse = ","), character(0))
  stretches_w_losses_num <- length(stretches_w_losses)
  stretches_w_losses_prop <- round(stretches_w_losses_num / num_stretches, 2)

  # Stretches with transitions
  stretches_w_transitions <- unique(c(stretches_w_losses, stretches_w_gains))
  stretches_w_transitions_str <- ifelse(length(stretches_w_transitions) > 0, paste0(stretches_w_transitions, collapse = ","), character(0))
  stretches_w_transitions_num <- length(stretches_w_transitions)
  stretches_w_transitions_prop <- round(stretches_w_transitions_num / num_stretches, 2)

  # Overall tarnsitions
  transitions <- c(gains, loss)
  transitions_str <-  ifelse(length(transitions) > 0, paste0(transitions, collapse = ","), character(0))
  transitions_num <- length(transitions)
  transitions_prop <- transitions_num / num_downstream_edges

  summary <- cbind.data.frame(num_stretches,
                             stretches_w_transitions = stretches_w_transitions_str, stretches_w_transitions_num,
                             stretches_w_gains = stretches_w_gains_str, stretches_w_gains_num, stretches_w_gains_prop,
                             stretches_w_losses = stretches_w_losses_str, stretches_w_losses_num, stretches_w_losses_prop,
                             num_downstream_edges,
                             transitions = transitions_str, transitions_num, transitions_prop,
                             gains = gains_str, gains_num, gains_prop,
                             loss = loss_str, loss_num, loss_prop)
  return(summary)
}
