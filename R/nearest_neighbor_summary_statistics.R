#' nearest_neighbor_summary_statistics: Generate summary statistics of nearest neighbor comparisons for isolate pairs
#'
#' This function will generate summary statistics for nearest neighbor comparisons
#'
#' @param nn_diff_df Dataframe with isolates and their nearest neighbors, generated from nearest_neighbor_analysis.
#' @param cat_vars Categorical variables. Must be in the comparison_df
#' @param cont_vars Continuous variables. Must be in the comparison_df
#' @param log_2 Whether to perform log-2 transformation on continuous variables
#' @return Dataframe with differences in the categorical and continuous variables between an isolate and their nearest neighbor
#' @export

nearest_neighbor_summary_statistics <- function (nn_diff_df, cat_vars, cont_vars, log_2){
  n <- nrow(nn_diff_df)
  # Category
  cat_summary <- function(var, nn_diff_df) {
    gain_prop <- sum(nn_diff_df[, var] == "gain")/n
    loss_prop <- sum(nn_diff_df[, var] == "loss")/n
    data <- cbind.data.frame(gain_prop, loss_prop) %>% `colnames<-`(paste0(var,
                                                                           c("_gain_prop", "_loss_prop")))
    return(data)
  }
  cat_vars_2 <- paste0(cat_vars, "_cat_diff")
  cat_sum_df <- lapply(cat_vars_2, FUN = function(x) {
    cat_summary(x, nn_diff_df)
  }) %>% do.call(cbind.data.frame, .)

  # Continuous
  cont_summary <- function(var, nn_diff_df) {
    string <- nn_diff_df[, var]
    median <- stats::median(string)
    mean <- mean(string)
    range <- range(string) %>% paste0(collapse = "-")
    increase_prop <- string %>% subset(. > 0) %>% length(.)/length(string)
    decrease_prop <- string %>% subset(. < 0) %>% length(.)/length(string)
    data <- cbind.data.frame(median, mean, range, increase_prop,
                             decrease_prop)
    colnames(data) <- paste0(var, "_", colnames(data))
    return(data)
  }
  if(log_2 == "yes"){
    cont_vars2_2 <- c(paste0(cont_vars, "_diff"),paste0(cont_vars,
                                                        "_log_2_diff"))
    cont_sum_df <- lapply(cont_vars2_2, FUN = function(x) {
      cont_summary(x, nn_diff_df)
    }) %>% do.call(cbind.data.frame, .)
  }
  else {
    cont_vars2_2 <- c(paste0(cont_vars, "_diff"))
    cont_sum_df <- lapply(cont_vars2_2, FUN = function(x) {
      cont_summary(x, nn_diff_df)
    }) %>% do.call(cbind.data.frame, .)
  }

  comps_sum <- cbind.data.frame(cat_sum_df, cont_sum_df)
  return(comps_sum)
}
