#' nearest_neighbor_analysis: Generate nearest neighbor comparisons for isolate pairs
#'
#' This function will compare the isolate to their nearest neighbor. Has functionality for categorical and continuous variables.
#'
#' This function will use the phylogenetic tree to identify an isolate's nearest neighbor.
#' @param comparison_df Dataframe with isolates and their nearest neighbors.
#' @param cat_vars Categorical variables. Must be in the comparison_df
#' @param cont_vars Continuous variables. Must be in the comparison_df
#' @param log_2 Whether to perform log-2 transformation on continuous variables
#' @return Dataframe with differences in the categorical and continuous variables between an isolate and their nearest neighbor
#' @export
nearest_neighbor_analysis <- function(comparison_df, cat_vars, cont_vars, log_2){
  nn_cat_diff <- function(df, var) {
    cat_diff <- ifelse(df[, var] == 1 & df[, paste0(var,
                                                    "_nn")] == 0, "gain", ifelse(df[, var] == 0 & df[,
                                                                                                     paste0(var, "_nn")] == 1, "loss", "same"))
    return(cat_diff)
  }
  cat_diff <- lapply(cat_vars, FUN = function(x) {
    nn_cat_diff(comparison_df, x)
  }) %>% do.call(cbind, .) %>% as.data.frame %>% `colnames<-`(paste0(cat_vars,
                                                                     "_cat_diff"))
  nn_cont_diff <- function(df, var, fc) {
    diff <- df[, var] - df[, paste0(var, "_nn")] %>% as.data.frame %>%
      `colnames<-`(paste0(var, "_diff"))
    if (log_2 == "yes"){
      log_2 <- log2(df[, var]) - log2(df[, paste0(var, "_nn")]) %>% as.data.frame %>%
        `colnames<-`(paste0(var, "_log_2_diff"))
      diff <- cbind.data.frame(diff, log_2)
    }
    else {
      diff <- diff
    }
    return(diff)
  }
  cont_diff <- lapply(cont_vars, FUN = function(x) {
    nn_cont_diff(comparison_df, x, fc = "yes")
  }) %>% do.call(cbind.data.frame, .)

  comps_df <- cbind.data.frame(isolate_no = comparison_df$isolate_no,nearest_neighbor=comparison_df$nearest_neighbor,cat_diff, cont_diff)
  return(comps_df)
}
