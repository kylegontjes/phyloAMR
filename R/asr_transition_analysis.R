#' asr_transition_analysis: Transition analysis
#'
#' This function calculates the transition statistics for an ancestral state reconstruction using their edge dataframe
#'
#' @param parent_child_df parent_child_df object from asr() object
#' @param node_states Whether this was done under joint or marginal reconstruction
#' @return A dataframe with data on the trait dynamics
#' @export
asr_transition_analysis <- function(parent_child_df,node_states="joint"){
  check_marginal_columns(node_states,parent_child_df)

  tip_data_df <- parent_child_df %>% subset(is.na(child_name) ==F)

  # Total data
  total_edges <- nrow(parent_child_df)

  # Transition data
  transitions <- parent_child_df[["transition"]] %>% sum
  losses <- parent_child_df[["loss"]] %>% sum
  losses_tip <- tip_data_df[["loss"]]%>% sum
  gains <- parent_child_df[["gain"]] %>% sum
  gains_tip <- tip_data_df[["gain"]] %>% sum

  # Continuation
  continuations <- parent_child_df[["continuation"]] %>% sum
  continuations_present <-  parent_child_df[["continuation_present"]]  %>% sum
  continuations_absent <-  parent_child_df[["continuation_absent"]]  %>% sum

  results <- cbind.data.frame(total_edges,transitions,gains,gains_tip,losses,losses_tip,continuations,continuations_present,continuations_absent)

  if(node_states =="marginal"){
    # Transition data
    transitions_high <- parent_child_df[["transition_high"]] %>% sum
    transitions_low <- parent_child_df[["transition_low"]] %>% sum
    losses_high <- parent_child_df[["loss_high"]] %>% sum
    losses_tip_high <- tip_data_df[["loss_high"]]%>% sum
    gains_high <- parent_child_df[["gain_high"]] %>% sum
    gains_tip_high <- tip_data_df[["gain_high"]] %>% sum

    # Continuation data
    continuations_unsure <- (parent_child_df[["continuation"]] == 1 & parent_child_df[["child_val"]]==0.5) %>% sum

    results <- cbind.data.frame(results,transitions_high,transitions_low,gains_high,gains_tip_high,losses_high,losses_tip_high,continuations_unsure) %>% select(total_edges,transitions,transitions_high,transitions_low,gains,gains_high,gains_tip,gains_tip_high,losses,losses_high,losses_tip,losses_tip_high,continuations,continuations_present,continuations_absent,continuations_unsure)
  }

  return(results)
}

check_marginal_columns <- function(node_states,parent_child_df){
  if(node_states == "marginal" & c(!"transition_high" %in% colnames(parent_child_df))){
    stop("Must have marginal columns of interest, ensure your parent_child_df is constructed accordingly")
  }
  if(node_states == "joint" & "transition_high" %in% colnames(parent_child_df)){
    stop("This might be a marginal reconstruction, specify node_state appropriately")
  }
}
