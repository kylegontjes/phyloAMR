#' Ancestral state reconstruction wrapper function
#'
#' Function to perform ancestral state reconstruction using corHMM. Output is parsed into a dataset with episodes of trait gain, loss, and continuation across the phylogeny
#'
#' @param df Dataframe with tip name (e.g., tip_name_variable) and phenotype/trait (e.g., trait) variables
#' @param tr Phylogeny of class "phylo"
#' @param tip_name_variable Name of variable containing tip names in df. Tip name variable must correspond to tip names in the tree.
#' @param trait Name of phenotype/trait variable in df
#' @param model Type of rate transition matrix. Options: equal rates model ("ER") or all rates different ("ARD"). An additional option, "MF", selects the best performing model using sample-size corrected Akaike information criterium (AICc). Default: ER
#' @param node_states Perform "joint" or "marginal" reconstruction. Default: joint
#' @param upper_bound Upper bound for likelihood search. Default: 1e50
#' @param lower_bound Lower bound for likelihood search. Default: 1e-9
#' @param confidence_threshold The confidence threshold to use when categorizing ancestral state inferences from marginal reconstruction. We highly suggest using a default of 0.5 (i.e., winner-takes-all), but permit modification to elevated values.
#' @return
#'   \describe{
#'     \item{corHMM_output}{corHMM output}
#'     \item{corHMM_model_statistics}{A dataframe containing the inferred rates, log-likelihood, AIC, and chosen transition model}
#'     \item{parent_child_df}{A dataframe of parent-child relationships with additional descriptive information}
#'     \item{node_states}{Text string indicating the chosen reconstruction method}
#'   }
#' @importFrom dplyr mutate
#' @export
asr <- function(df, tr, tip_name_variable, trait, model = "ER", node_states = "joint", upper_bound = 1e50, lower_bound = 1e-9, confidence_threshold = NULL) {
  # Check if phenotype is 0,1
  check_trait(df[[trait]])
  # Check if tree is rooted
  check_tree(tr)
  # Check if phenotype and tree length is the same
  check_trait_tree_length(trait_var = df[[trait]], df_tips = df[[tip_name_variable]], tree = tr)
  # Check if states appropriately chosen
  check_joint_confidence_value(node_states, confidence_threshold)

  # Order dataframe properly, if not the predictions will be wrong
  df <- df[match(tr$tip.label, df[[tip_name_variable]]), ]

  # Run corHMM to estimate hidden rates
  ## Check if you want to use the best model
  if (model == "MF") {
    cat("Performing model finding to use best fitting model \n")
    MF_results <- find_best_asr_model(df = df, tr = tr, tip_name_variable = tip_name_variable, trait = trait, node_states = node_states, upper_bound = upper_bound, lower_bound = lower_bound)
    corHMM_output <- MF_results[["corHMM_output_best_model"]]
  } else {
    corHMM_output <- corHMM::corHMM(phy = tr, data = df[, c(tip_name_variable, trait)], rate.cat = 1, model = model, node.states = node_states, upper.bound = upper_bound, lower.bound = lower_bound)
  }

  # Check if rates are at min or max bound
  check_rates_at_local_max(corHMM_output, upper_bound = upper_bound, lower_bound = lower_bound)

  # Model stats
  if (model == "MF") {
    corHMM_model_statistics <- MF_results[["model_statistics"]]
    # Remove model finder object
    rm(MF_results)
  } else {
    corHMM_model_statistics <- characterize_asr_model(corHMM_output)
  }

  # Get parent child data using the edge matrix, ancestral state predictions, and trait data
  outcome_str <- df[, trait]
  names(outcome_str) <- df[[tip_name_variable]]
  parent_child_df <- get_parent_child_data(tr = tr, ancestral_states = corHMM_output$states, trait_data = outcome_str, node_states = node_states, confidence_threshold = confidence_threshold)

  # Annotate parent child data with transition data (i.e., gain, loss, and continuation)
  parent_child_df <- get_continuation_data(parent_child_df = parent_child_df, node_states = node_states)

  # Results object
  asr_output  <- list(corHMM_output = corHMM_output,
                      corHMM_model_statistics = corHMM_model_statistics,
                      parent_child_df = parent_child_df,
                      node_states = node_states)

  return(asr_output)
}

#' Get parent child data from ancestral state reconstruction
#'
#' Function to generate an annotated edge matrix with ancestral and tip states
#'
#' @param tr Phylogenetic tree of class "phylo"
#' @param ancestral_states Ancestral states inferred by corHMM
#' @param trait_data Named phenotype/trait string
#' @param node_states Perform "joint" or "marginal" reconstruction. Default: joint
#' @param confidence_threshold The confidence threshold to use for marginal state reconstruction. Suggested value: 0.875.
#' @return
#'   \describe{
#'     \item{edge}{Edge matrix with ancestral and tip states}
#'   }
#' @export
get_parent_child_data <- function(tr, ancestral_states, trait_data, node_states, confidence_threshold = NULL) {
  # Tree edge info
  edge <- tr$edge %>% as.data.frame
  colnames(edge) <- c("parent", "child")

  # Prediction data
  if (node_states == "marginal") {
    ancestral_states <- as.data.frame(ancestral_states) %>% `colnames<-`(levels(as.factor(trait_data)))
    ancestral_states[, "pred"] <- ifelse(ancestral_states[, 2] > confidence_threshold, 1, ifelse(ancestral_states[, 1] > confidence_threshold, 0, 0.5))
  } else if (node_states == "joint") {
    state_coding <- levels(as.factor(trait_data))
    names(state_coding) <- seq_along(state_coding)
    ancestral_states <- data.frame(pred = ancestral_states)
    ancestral_states$pred <- dplyr::recode(ancestral_states$pred, !!!as.list(state_coding)) %>% as.numeric
  } else {
    stop("This tool is only optimized for marginal and joint corHMM")
  }

  # Add node data
  ancestral_states[, "node"] <-  seq_len(nrow(ancestral_states)) + (nrow(ancestral_states) + 1)
  num_edges <- nrow(edge)
  internal_nodes <- c(length(trait_data) + 1):(unlist(edge) %>% max)

  # Dataframe construction
  edge$parent_val <- rep(NA, num_edges)
  edge$child_val <- rep(NA, num_edges)

  # Prediction values for the parent & child
  for (i in internal_nodes) {
    edge[edge$parent == i, "parent_val"] <- ancestral_states[ancestral_states$node == i, "pred"]
    edge[edge$child == i, "child_val"] <- ancestral_states[ancestral_states$node == i, "pred"]
  }

  # Internal nodes
  if (node_states == "marginal") {
    edge$child_mle <- rep(NA, num_edges)
    edge$parent_mle <- rep(NA, num_edges)
    for (i in internal_nodes){
      edge[edge$parent == i, "parent_mle"] <- ancestral_states[ancestral_states$node == i, levels(as.factor(trait_data))[2]] * 100
      edge[edge$child == i, "child_mle"] <- ancestral_states[ancestral_states$node == i, levels(as.factor(trait_data))[2]] * 100
    }
  }

  # Values and name for the tip
  for (i in 1:(min(internal_nodes) - 1)) {
    edge[edge$child == i, "child_name"] <- tr$tip.label[[i]]
    edge[edge$child == i, "child_val"] <- subset(trait_data, names(trait_data) == edge[edge$child == i, "child_name"])
  }

  return(edge)
}

#' Annotate parent child data from ancestral state reconstructionwith transition data
#'
#' Function to generate an annotated edge matrix with transition data
#'
#' @param parent_child_df Edge matrix with ancestral states
#' @param node_states Whether to perform "joint" or "marginal" reconstruction. Default: joint
#' @return
#'   \describe{
#'     \item{parent_child_df}{Final parent child dataset with transition data}
#'   }
#' @export
get_continuation_data <- function(parent_child_df, node_states) {
  if (node_states == "joint") {
    parent_child_df <- parent_child_df %>% dplyr::mutate(transition = ifelse(parent_val != child_val, 1, 0),
                                                         gain = ifelse(child_val == 1 & parent_val == 0, 1, 0),
                                                         loss = ifelse(child_val == 0 & parent_val == 1, 1, 0),
                                                         continuation  =  ifelse(parent_val == child_val, 1, 0),
                                                         continuation_present = ifelse(continuation == 1 & child_val == 1, 1, 0),
                                                         continuation_absent = ifelse(continuation == 1 & child_val == 0, 1, 0))
  } else if (node_states == "marginal") {
    parent_child_df <- parent_child_df %>% dplyr::mutate(transition = ifelse(parent_val != child_val, 1, 0),
                                                         transition_high = ifelse(parent_val == 0 & child_val == 1 | parent_val == 1 & child_val == 0, 1, 0),
                                                         transition_low = ifelse(parent_val == 0.5 & child_val == 1 | parent_val == 0.5 & child_val == 0, 1, 0),
                                                         gain = ifelse(child_val == 1 & parent_val == 0.5 | child_val == 1 & parent_val == 0, 1, 0),
                                                         gain_high  =  ifelse(child_val == 1 & parent_val == 0, 1, 0),
                                                         gain_low = ifelse(child_val == 1 & parent_val == 0.5, 1, 0),
                                                         loss = ifelse(child_val == 0 & parent_val == 1 | child_val == 0 & parent_val == 0.5, 1, 0),
                                                         loss_high  =  ifelse(child_val == 0 & parent_val == 1, 1, 0),
                                                         loss_low  =  ifelse(child_val == 0 & parent_val == 0.5, 1, 0),
                                                         continuation  =  ifelse(parent_val == child_val, 1, 0),
                                                         continuation_high  =  ifelse(parent_val == 1 & child_val == 1 | parent_val == 0 & child_val == 0, 1, 0),
                                                         continuation_low =  ifelse(parent_val == 0.5 & child_val == 0.5, 1, 0),
                                                         continuation_present = ifelse(continuation_high == 1 & child_val == 1, 1, 0),
                                                         continuation_absent = ifelse(continuation_high == 1 & child_val == 0, 1, 0))
  }
  return(parent_child_df)
}
