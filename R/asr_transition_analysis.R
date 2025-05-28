#' Enumeration of trait transitions using ancestral state reconstruction
#'
#' Calculate transition statistics (e.g., gain, loss, and continuation) of an ancestral state reconstruction model.
#'
#' @param parent_child_df parent_child_df object returned sfrom asr() function
#' @param node_states If 'joint' or 'marginal' reconstruction was used
#' @return Dataframe with phylogenetic transition statistics for a trait
#' @export
asr_transition_analysis <- function(parent_child_df, node_states = "joint") {
  # Checks
  ## Does appropriate columns exist when requesting statistics for marginal reconstruction
  check_marginal_columns(node_states, parent_child_df)

  # Subset parent child dataframe to just include tips
  tip_data_df <- subset(parent_child_df, !is.na(parent_child_df$child_name))

  # Total data
  total_edges <- nrow(parent_child_df)

  # Transition data
  transitions <- sum(parent_child_df[["transition"]])
  losses <- sum(parent_child_df[["loss"]])
  losses_tip <- sum(tip_data_df[["loss"]])
  gains <- sum(parent_child_df[["gain"]])
  gains_tip <- sum(tip_data_df[["gain"]])

  # Continuation
  continuations <- sum(parent_child_df[["continuation"]])
  continuations_present <- sum(parent_child_df[["continuation_present"]])
  continuations_absent <-  sum(parent_child_df[["continuation_absent"]])

  # Parent data for summary stats
  parents_w_trait <- sum(parent_child_df[["parent_value"]])
  parents_wo_trait <- total_edges - parents_w_trait

  # Summary Statistics
  gain_frequency <- round(gains / parents_wo_trait * 100, 2)
  loss_frequency <- round(losses / parents_w_trait * 100, 2)
  continuation_present_frequency <- round(continuations_present / total_edges * 100, 2)
  continuation_absent_frequency <- round(continuations_absent / total_edges * 100, 2)

  results <- cbind.data.frame(total_edges, transitions, gains, gains_tip, losses, losses_tip,
                              continuations, continuations_present, continuations_absent,
                              gain_frequency, loss_frequency,
                              continuation_present_frequency, continuation_absent_frequency)

  if (node_states == "marginal") {
    # Transition data
    transitions_high <- sum(parent_child_df[["transition_high"]])
    transitions_low <- sum(parent_child_df[["transition_low"]])
    losses_high <- sum(parent_child_df[["loss_high"]])
    losses_tip_high <- sum(tip_data_df[["loss_high"]])
    gains_high <- sum(parent_child_df[["gain_high"]])
    gains_tip_high <- sum(tip_data_df[["gain_high"]])

    # Continuation data
    continuations_unsure <- sum(parent_child_df[["continuation"]] == 1 & parent_child_df[["child_value"]] == 0.5)

    # Compile all information in a results dataframe
    results <- cbind.data.frame(results, transitions_high, transitions_low, gains_high, gains_tip_high, losses_high, losses_tip_high, continuations_unsure)
    results <- results  %>% select(total_edges, transitions, transitions_high, transitions_low,
                                   gains, gains_high, gains_tip, gains_tip_high,
                                   losses, losses_high, losses_tip, losses_tip_high,
                                   continuations, continuations_present, continuations_absent, continuations_unsure,
                                   gain_frequency, loss_frequency, continuation_present_frequency, continuation_absent_frequency)
  }

  return(results)
}
