#' Analysis of synchronous transitions of two traits
#'
#' Synchronous gain and loss of two traits at a spot on tree
#'
#' @param comparitor_parent_child_df Parent child dataset for a comparitor trait, such as a genotype or a different trait/phenotype
#' @param trait_parent_child_df Parent child dataset for a trait of interest, such as a trait/phenotype
#' @param node_states Joint or marginal reconstruction
#' @param confidence Whether to use high or low confidence transition nodes when node_states are marginal
#' @return Synchronous gain and loss events of two traits
#' @export
synchronous_transitions <- function(comparitor_parent_child_df, trait_parent_child_df, node_states = "joint", confidence = NULL) {
  # Get gain and loss nodes for the comparitor and trait
  if (node_states == "joint") {
    trait_gains <- trait_parent_child_df[trait_parent_child_df$gain == 1, "child"]
    comparitor_gains <- comparitor_parent_child_df[comparitor_parent_child_df$gain == 1, "child"]
    trait_losses <- trait_parent_child_df[trait_parent_child_df$loss == 1, "child"]
    comparitor_losses <- comparitor_parent_child_df[comparitor_parent_child_df$loss == 1, "child"]
  } else {
    # Get variable to consider
    gain_variable <- ifelse(confidence == "high", "gain_high", "gain")
    loss_variable <- ifelse(confidence == "high", "loss_high", "loss")
    # Identify
    trait_gains <- trait_parent_child_df[trait_parent_child_df[[gain_variable]] == 1, "child"]
    comparitor_gains <- comparitor_parent_child_df[comparitor_parent_child_df[[gain_variable]] == 1, "child"]
    trait_losses <- trait_parent_child_df[trait_parent_child_df[[loss_variable]] == 1, "child"]
    comparitor_losses <- trait_parent_child_df[trait_parent_child_df[[loss_variable]] == 1, "child"]
  }

  # Trait
  num_trait_gains <- length(trait_gains)
  num_trait_losses <- length(trait_losses)

  # Synchronous gains
  synchronous_gains  <- trait_gains[which(trait_gains %in% comparitor_gains)]
  synchronous_gains <- sort(as.numeric(synchronous_gains))
  synchronous_gains_str <-  if (length(synchronous_gains) == 0) "" else paste0(synchronous_gains, collapse = ",")
  synchronous_gains_num <- length(synchronous_gains)
  synchronous_gains_prop <- round(synchronous_gains_num / num_trait_gains, 2)

  # Synchronous gains and losses
  synchronous_gain_loss <- trait_gains[which(trait_gains %in%  comparitor_losses)]
  synchronous_gain_loss <- sort(as.numeric(synchronous_gain_loss))
  synchronous_gain_loss_str <- if (length(synchronous_gain_loss) == 0) "" else paste0(synchronous_gain_loss, collapse = ",")
  synchronous_gain_loss_num <- length(synchronous_gain_loss)
  synchronous_gain_loss_prop <- round(synchronous_gain_loss_num / num_trait_gains, 2)

  # Synchronous losses
  synchronous_losses <- trait_losses[which(trait_losses %in%  comparitor_losses)]
  synchronous_losses <- sort(as.numeric(synchronous_losses))
  synchronous_losses_str <- if (length(synchronous_losses) == 0) "" else paste0(synchronous_losses, collapse = ",")
  synchronous_losses_num <- length(synchronous_losses)
  synchronous_losses_prop <- round(synchronous_losses_num / num_trait_losses, 2)

  # Synchronous loss gain
  synchronous_loss_gain <- trait_losses[which(trait_losses %in%  comparitor_gains)]
  synchronous_loss_gain <- sort(as.numeric(synchronous_loss_gain))
  synchronous_loss_gain_str <- if (length(synchronous_loss_gain) == 0) "" else paste0(synchronous_loss_gain, collapse = ",")
  synchronous_loss_gain_num <- length(synchronous_loss_gain)
  synchronous_loss_gain_prop <- round(synchronous_loss_gain_num / num_trait_losses, 2)

  # Number synchronous transitions
  synchronous_transitions_num <- sum(synchronous_gains_num, synchronous_gain_loss_num, synchronous_losses_num, synchronous_loss_gain_num)

  summary <- data.frame(num_trait_gains, synchronous_transitions_num, synchronous_gains = synchronous_gains_str, synchronous_gains_num, synchronous_gains_prop,
                        synchronous_gain_loss = synchronous_gain_loss_str, synchronous_gain_loss_num, synchronous_gain_loss_prop,
                        num_trait_losses, synchronous_losses = synchronous_losses_str, synchronous_losses_num, synchronous_losses_prop,
                        synchronous_loss_gain = synchronous_loss_gain_str, synchronous_loss_gain_num, synchronous_loss_gain_prop, stringsAsFactors = TRUE)
  return(summary)
}
