#' asr_transition_analysis: Transition analysis
#'
#' This function calculates the transition statistics for ancestral state
#'
#' @param parent_child_df parent_child_df object from asr() object
#' @param node_states Whether this was done under joint or marginal reconstruction
#' @return A dataframe with data on the trait dynamics
#' @export
asr_transition_analysis <- function(parent_child_df,node_states="joint"){
  tip_data_df <- parent_child_df %>% subset(is.na(child_name) ==F)

  # Total data
  total_edges <- nrow(parent_child_df)

  # Transition data
  transitions <- parent_child_df[["transition_any"]] %>% sum
  losses <- parent_child_df[["loss_any"]] %>% sum
  losses_tip <- tip_data_df[["loss_any"]]%>% sum
  gains <- parent_child_df[["gain_any"]] %>% sum
  gains_tip <- tip_data_df[["gain_any"]] %>% sum

  results <- cbind.data.frame(total_edges,transitions,gains,gains_tip,losses,losses_tip)

  if(node_states =="marginal"){
    transitions_high <- parent_child_df[["transition_high"]] %>% sum
    losses_high <- parent_child_df[["loss_high"]] %>% sum
    losses_tip_high <- tip_data_df[["loss_high"]]%>% sum
    gains_high <- parent_child_df[["gain_high"]] %>% sum
    gains_tip_high <- tip_data_df[["gain_high"]] %>% sum
    results <- cbind.data.frame(results,transitions_high,gains_high,gains_tip_high,losses_high,losses_tip_high) %>% select(transitions,transitions_high,gains,gains_high,gains_tip,gains_tip_high,losses,losses_high,losses_tip,losses_tip_high)
  }

  return(results)
}

