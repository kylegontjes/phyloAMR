#' nearest_neighbor_summary_statistics: Generate summary statistics of nearest neighbor comparisons for isolate pairs
#'
#' This function will generate summary statistics for nearest neighbor comparisons
#'
#' @param nn_diff_df Dataframe with isolates and their nearest neighbors, generated from nearest_neighbor_analysis.
#' @param categorical_vars Categorical variables. Must be in the comparison_df
#' @param binary_vars Binary variables. Must be in the comparison_df AND coded 0 and 1.
#' @param continuous_vars Continuous variables. Must be in the comparison_df
#' @param continuous_vars Continuous variables. Must be in the comparison_df
#' @param log_2 Boolean (TRUE/FALSE) Whether to perform log-2 transformation on continuous variables
#' @return Dataframe with differences in the categorical and continuous variables between an isolate and their nearest neighbor
#' @export

nearest_neighbor_summary_statistics <- function(nn_diff_df, categorical_vars = NULL, binary_vars = NULL, continuous_vars = NULL, log_2 = FALSE) {
  n <- nrow(nn_diff_df)
  # Create summary dataset
  comps_sum <- cbind.data.frame(n = n)

  # Categorical differences
  if (is.null(categorical_vars) == FALSE) {
    continuous_vars_names <- paste0(categorical_vars, "_cat_diff")
    cat_sum_df <- lapply(continuous_vars_names, FUN = function(x) {
      cat_summary(x, nn_diff_df = nn_diff_df, n = n)
    }) %>% do.call(cbind.data.frame, .)
    comps_sum <- cbind(comps_sum, cat_sum_df)
  }

  # Binary differences
  if (is.null(binary_vars) == FALSE) {
    binary_vars_names <- paste0(binary_vars, "_binary_diff")
    binary_sum_df <- lapply(binary_vars_names, FUN = function(x) {
      binary_summary(x, nn_diff_df = nn_diff_df, n = n)
    }) %>% do.call(cbind.data.frame, .)
    comps_sum <- cbind(comps_sum, binary_sum_df)
  }

  # Continuous
  if (is.null(continuous_vars) == FALSE) {
    if (log_2 == TRUE) {
      continuous_vars_names <- c(paste0(continuous_vars, "_diff"), paste0(continuous_vars, "_log_2_diff"))
      cont_sum_df <- lapply(continuous_vars_names, FUN = function(x) {
        cont_summary(x, nn_diff_df = nn_diff_df, n = n)
      }) %>% do.call(cbind.data.frame, .)
    } else {
      continuous_vars_names <- c(paste0(continuous_vars, "_diff"))
      cont_sum_df <- lapply(continuous_vars_names, FUN = function(x) {
        cont_summary(x, nn_diff_df = nn_diff_df, n = n)
      }) %>% do.call(cbind.data.frame, .)
  }
    comps_sum <- cbind(comps_sum, cont_sum_df)
  }

  return(comps_sum)
}

cat_summary <- function(var, nn_diff_df, n) {
  difference_prop <- sum(nn_diff_df[, var] != "same") / n
  data <- cbind.data.frame(difference_prop) %>% `colnames<-`(paste0(var,  c("_difference_prop")))
  return(data)
}

binary_summary <- function(var, nn_diff_df, n) {
  gain_prop <- sum(nn_diff_df[, var] == "gain") / n
  loss_prop <- sum(nn_diff_df[, var] == "loss") / n
  data <- cbind.data.frame(gain_prop, loss_prop) %>% `colnames<-`(paste0(var,  c("_gain_prop", "_loss_prop")))
  return(data)
}

cont_summary <- function(var, nn_diff_df, n) {
  string <- nn_diff_df[, var]
  median <- stats::median(string)
  mean <- mean(string)
  range <- range(string) %>% paste0(collapse = "-")
  increase_prop <- string %>% subset(. > 0) %>% length(.) / length(string)
  decrease_prop <- string %>% subset(. < 0) %>% length(.) / length(string)
  data <- cbind.data.frame(median, mean, range, increase_prop, decrease_prop)
  colnames(data) <- paste0(var, "_", colnames(data))
  return(data)
}
