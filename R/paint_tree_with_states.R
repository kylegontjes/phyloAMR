#' Paint phylogenetic tree with ancestral and tip states
#'
#' Function to paint a tree with ancestral and tip statess
#'
#' @param parent_child_df Parent child dataframe from asr()
#' @param tr Phylogenetic tree
#' @param scale_breaks Values for breaks (marginal = c(1,0.5,0); joint = c(1,0))
#' @param scale_labels Labels for present and absent, marginal should be c('Present','Unsure','Absent')
#' @param scale_colors Colors for states
#' @param legend_name Name for legend
#' @return Tree with painted states
#' @export
paint_tree_with_states <- function(parent_child_df, tr, scale_breaks = c(1, 0), scale_labels = c("Present", "Absent"), scale_colors = c("red", "black"), legend_name = "Trait") {
  # Create scale
  edge_color <-  ggtree::scale_color_manual(breaks = scale_breaks, values = scale_colors, labels = scale_labels, name = legend_name)

  # Get number of nodes
  num_nodes <- nrow(tr$edge) + 1
  root_node <- length(tr$tip.label) + 1

  # Get root state
  root_val <- unique(parent_child_df[parent_child_df$parent == root_node, "parent_val"])

  # Get outcome and parent strings
  ## Order parent child dataframe
  parent_child_df <- parent_child_df[order(parent_child_df$child, decreasing = FALSE), ]
  number_tips <- length(tr$tip.label)
  ## Get outcome and parent string
  outcome_str <- parent_child_df[parent_child_df$child <= number_tips, "child_val"]
  parent_str <- parent_child_df[parent_child_df$child > number_tips, "child_val"]

  # Merge as dataframe
  state_df <- data.frame(node = seq_len(num_nodes), state = c(outcome_str, root_val, parent_str))
  tr <- dplyr::full_join(tr,state_df)

  # Tree
  painted_tree <- ggtree::ggtree(tr, aes(color = as.factor(state))) + edge_color
  return(painted_tree)
}
