#' nearest_neighbor_analysis: Generate nearest neighbor comparisons for isolate pairs
#'
#' This function will compare the isolate to their nearest neighbor generated from phyloAMR::nearest_neighbor_algorithm(). Has functionality for categorical and continuous variables.
#'
#' This function will use the phylogenetic tree to identify an isolate's nearest neighbor.
#' @param nearest_neighbor_df Dataframe with isolates and their nearest neighbors. Can work on
#' @param categorical_vars Categorical variables. Must be in the nearest_neighbor_df
#' @param binary_vars Binary variables. Must be in the nearest_neighbor_df AND coded 0 and 1.
#' @param continuous_vars Continuous variables. Must be in the nearest_neighbor_df
#' @param log_2 Boolean (TRUE/FALSE) Whether to perform log-2 transformation on continuous variables
#' @return Dataframe with differences in the categorical and continuous variables between an isolate and their nearest neighbor
#' @export
nearest_neighbor_analysis <- function(nearest_neighbor_df, categorical_vars = NULL, binary_vars = NULL, continuous_vars = NULL, log_2 = FALSE) {
  # Create comparison dataframe
  comps_df <- nearest_neighbor_df[, c("tip_name","nearest_neighbor","comparison_feature","isolate_comparison_value","nearest_neighbor_comparison_value")]

  # Categorical differences
  if (is.null(categorical_vars) == FALSE) {
    cat_diff <- lapply(categorical_vars, FUN = function(x) {
      nn_cat_diff(nearest_neighbor_df, x)
    })
    cat_diff <- do.call(cbind.data.frame, cat_diff)
    colnames(cat_diff) <- paste0(categorical_vars, "_categorical_diff")
    comps_df <- cbind(comps_df, cat_diff)
  }

  # Binary differences
  if (is.null(binary_vars) == FALSE) {
    binary_diff <- lapply(binary_vars, FUN = function(x) {
      nn_binary_diff(nearest_neighbor_df, x)
    })
    binary_diff <- do.call(cbind.data.frame, binary_diff)
    colnames(binary_diff) <- paste0(binary_vars, "_binary_diff")
    comps_df <- cbind(comps_df, binary_diff)
  }

  # Continuous differences
  if (is.null(continuous_vars) == FALSE) {
    cont_diff <- lapply(continuous_vars, FUN = function(x) {
      nn_cont_diff(nearest_neighbor_df, x, log_2 = log_2)
    })
    cont_diff <- do.call(cbind.data.frame, cont_diff)
    comps_df <- cbind(comps_df, cont_diff)
  }

  return(comps_df)
}

# Categorical differences
nn_cat_diff <- function(df, var) {
  cat_diff <- ifelse(df[, var] ==  df[, paste0(var, "_nn")], "same", "different")
  return(cat_diff)
}

# Binary differences
nn_binary_diff <-  function(df, var) {
  bin_diff <- ifelse(df[, var] == 1 & df[, paste0(var, "_nn")] == 0, "gain", ifelse(df[, var] == 0 & df[, paste0(var, "_nn")] == 1, "loss", "same"))
  return(bin_diff)
}

# Continuous differences
nn_cont_diff <- function(df, var, log_2) {
  diff <- df[, var] - df[, paste0(var, "_nn")]
  diff <- as.data.frame(diff)
  colnames(diff) <- paste0(var, "_diff")

  if (log_2 == TRUE) {
    log_2_diff <- log2(df[, var]) - log2(df[, paste0(var, "_nn")])
    log_2_diff <- as.data.frame(log_2_diff)
    colnames(log_2_diff) <- paste0(var, "_log_2_diff")
    diff <- cbind.data.frame(diff, log_2_diff)
  } else {
    diff <- diff
  }
  return(diff)
}
