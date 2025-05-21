#' downstream_gain_loss: Analysis of downstream gain and loss of comparitor traits
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
downstream_gain_loss <- function(comparitor_parent_child_df, trait_parent_child_df, tr, node_states = "joint", confidence = NULL) {
  stretches <- get_trait_traces_on_tree(parent_child_df = trait_parent_child_df, tr = tr, node_states = node_states)
  downstream_changes <- get_gain_loss_on_stretches(comparitor_parent_child_df = comparitor_parent_child_df, stretches = stretches, node_states = node_states, confidence = confidence)
  return(downstream_changes)
}

get_trait_traces_on_tree <- function(parent_child_df, tr, node_states) {
  all_possible_paths <- ape::nodepath(tr)

  # Get gain location
  if (node_states == "joint") {
    gains <- parent_child_df %>% subset(gain == 1) %>% .$child
  } else {
    gains <- parent_child_df %>% subset(get(paste0("gain_", confidence)) == 1) %>% .$child
  }

  # Decipher tip and ancestral gains
  tip_gains <- subset(gains, gains < min(tr$edge[1, ]))
  ancestral_gains <- subset(gains, !gains %in% tip_gains)

  gain_paths <- lapply(all_possible_paths, FUN = function(x, gains) {
    didit <- x %in% gains
    if (sum(didit) == 0) {
      string <- NA
    } else {
      didit_pos <- which(didit)
      string <- x[didit_pos:length(x)]
    }
    return(string)
  }, gains = ancestral_gains) %>% subset(is.na(.) == FALSE)

  trait_found <- unique(c(subset(parent_child_df, child_val == 1) %>% .$child, subset(parent_child_df, parent_val == 1) %>% .$parent))

  paths_w_trait <- lapply(gain_paths, FUN = function(x, trait_found) {
    results <- x %in% trait_found
    names(results) <- x
    return(results)
  }, trait_found = trait_found)

  find_true_stretches <- function(x) {
    rle_result <- rle(x)  # Run-length encoding
    nodes <- names(x)
    starts <- cumsum(c(1, head(rle_result$lengths, -1)))  # Start positions
    stretches <- data.frame(startnode = names(starts[rle_result$values]), start = starts[rle_result$values],
                            end = starts[rle_result$values] + rle_result$lengths[rle_result$values] - 1) %>% subset(start == 1)

    start <- stretches[1, "start"]
    end <- stretches[1, "end"]
    paths <- as.vector(nodes[start:end])
    return(paths)
  }
  true_stretches <- lapply(paths_w_trait, find_true_stretches) %>% unique
  return(true_stretches)
}

get_gain_loss_on_stretches <- function(comparitor_parent_child_df, stretches, node_states, confidence) {
  if (node_states == "joint") {
    comparitor_gains <- comparitor_parent_child_df  %>% subset(gain == 1) %>% .$child
    comparitor_losses <- comparitor_parent_child_df  %>% subset(loss == 1) %>% .$child
  } else {
    comparitor_gains <- comparitor_parent_child_df  %>% subset(get(paste0("gain_", confidence)) == 1) %>% .$child
    comparitor_losses <- comparitor_parent_child_df  %>% subset(get(paste0("loss_", confidence)) == 1) %>% .$child
  }

  get_downstream_nodes <- function(stretches) {
    lapply(stretches, FUN = function(x) {
      x[-1]
      })
  }

  downstream <- get_downstream_nodes(stretches)
  downstream_positions <- stretches %>% unlist %>% unique %>% sort

  num_downstream_edges <- length(downstream_positions)
  gains <- downstream_positions[which(downstream_positions %in% comparitor_gains)] %>% as.numeric %>% sort
  gains_str <- paste0(gains, collapse = ",")
  gains_num <- length(gains)
  gains_prop <- gains_num / num_downstream_edges
  loss <- downstream_positions[which(downstream_positions %in% comparitor_losses)] %>% as.numeric %>% sort
  loss_str <- paste0(loss, collapse = ",")
  loss_num <- length(loss)
  loss_prop <- loss_num / num_downstream_edges

  # Collated lists
  get_gain_paths <- function(stretches) {
    gains <- lapply(stretches, FUN = function(x) {x[1]}) %>% unlist %>% unique %>% sort
    # group by stretches with same start
    paths_from_gain <- lapply(gains, FUN = function(gain) {
      stretches[lapply(stretches, FUN = function(y) {gain %in% y}) %>% unlist %>% which] %>% unlist %>% as.numeric %>% sort(decreasing = TRUE) %>% unique()
    })

    downstream_merged_paths <- lapply(paths_from_gain, FUN = function(x) {
      without_downstream_but_merged <- subset(x, !x %in% gains)
      return(without_downstream_but_merged)})
    names(downstream_merged_paths) <- lapply(paths_from_gain, FUN = function(x) {
      subset(x, x %in% gains)
    }) %>% unlist
    return(downstream_merged_paths)
  }
  merged_paths <- get_gain_paths(stretches = stretches)
  # Num paths
  num_stretches <- length(merged_paths)
  #Gains
  stretches_w_gains <- if(length(merged_paths) == 0){
    NULL
    } else {
    lapply(merged_paths, FUN = function(x) {
    x[which(x %in% comparitor_gains)] %>% sum %>% {ifelse(. > 0, TRUE, FALSE)}
  }) %>% unlist %>% which %>% names %>% as.numeric %>% sort
  }
  stretches_w_gains_str <- ifelse(length(stretches_w_gains) > 0, paste0(stretches_w_gains, collapse = ","), "")
  stretches_w_gains_num <- length(stretches_w_gains)
  stretches_w_gains_prop <- stretches_w_gains_num / num_stretches
  # Losses
  stretches_w_losses <- if(length(merged_paths) == 0){
    NULL
  } else {
    lapply(merged_paths, FUN = function(x) {
    x[which(x %in% comparitor_losses)] %>% sum %>% {ifelse(. > 0, TRUE, FALSE)}
  }) %>% unlist %>% which %>% names %>% as.numeric %>% sort
  }
  stretches_w_losses_str <- ifelse(length(stretches_w_losses) > 0, paste0(stretches_w_losses, collapse = ","), "")
  stretches_w_losses_num <- length(stretches_w_losses)
  stretches_w_losses_prop <- stretches_w_losses_num / num_stretches

  stretches_w_transitions <- unique(c(stretches_w_losses,stretches_w_gains))
  stretches_w_transitions_str <- ifelse(length(stretches_w_transitions) > 0, paste0(stretches_w_transitions, collapse = ","), "")
  stretches_w_transitions_num <- length(stretches_w_transitions %>% subset(.!=''))
  stretches_w_transitions_prop <- stretches_w_transitions_num / num_stretches

  transitions <- c(gains,loss)
  transitions_str <-  ifelse(length(transitions) > 0, paste0(transitions, collapse = ","), "")
  transitions_num <- length(transitions)
  transitions_prop <- transitions_num / num_downstream_edges

  summary <- cbind.data.frame(num_stretches,
                             stretches_w_transitions = stretches_w_transitions_str, stretches_w_transitions_num,
                             stretches_w_gains = stretches_w_gains_str, stretches_w_gains_num, stretches_w_gains_prop,
                             stretches_w_losses = stretches_w_losses_str, stretches_w_losses_num, stretches_w_losses_prop,
                             num_downstream_edges,
                             transitions = transitions_str, transitions_num,transitions_prop,
                             gains = gains_str, gains_num, gains_prop,
                             loss = loss_str, loss_num, loss_prop)
  return(summary)
}
