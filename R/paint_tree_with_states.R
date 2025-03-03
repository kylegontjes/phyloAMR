#' paint_tree_with_states
#'
#' Function to paint a tree with anestral states
#'
#' @param parent_child_df Parent child dataframe from asr()
#' @param tr Phylogenetic tree
#' @param scale_labels Labels for present and absent
#' @param legend_name Name for legend
#' @return Tree with painted states
#' @export
paint_tree_with_states <- function(parent_child_df,tr,scale_breaks=c(1,0),scale_labels = c("Present","Absent"),scale_colors=c("red","black"),legend_name = "Trait"){
  edge_color <-  scale_color_manual(breaks=scale_breaks,values=scale_colors,labels=scale_labels,name=legend_name)
  # Get outcome and parent string
  outcome_str = parent_child_df %>% subset(.,child <= length(tr$tip.label)) %>% arrange(child)  %>% .$child_val
  parent_str = parent_child_df %>% subset(.,child > length(tr$tip.label)) %>% arrange(child) %>% .$child_val
  # Get data
  d = data.frame(node=seq_len(nrow(parent_child_df)),state=c(outcome_str,parent_str))
  # Tree
  painted_tree <- ggtree(tr) %<+% d + aes(color=as.factor(state)) + edge_color
  return(painted_tree)
}
