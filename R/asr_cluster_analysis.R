#' asr_cluster_analysis: Cluster analysis
#'
#' This function calculates the clustering statistics (e.g., number of singletons, clusters, and no features)
#'
#' @param tip_data_df Tip data frame object from asr_cluster_detection()
#' @return A dataframe with data on the phylogenetics of a trait
#' @export
asr_cluster_analysis <- function(tip_data_df) {
  clustering <- tip_data_df[["asr_cluster"]]

  sum_clustering <-  table(clustering)
  num_isolates  <-  length(clustering)

  # Characterize present
  present_isolates <- sum_clustering[grepl("cluster|singleton", names(sum_clustering))]
  present <- sum(present_isolates)
  singletons <- sum(present_isolates[grepl("singleton",names(present_isolates))], sum(grepl("1pt",names(present_isolates))))
  singleton_isolates <- sum(present_isolates[grepl("singleton|1pt",names(present_isolates))])
  if (singleton_isolates == present) {
    clusters <- 0
    cluster_isolates <-  0
    cluster_size_median <- 0
    cluster_size_range <- NA
    cluster_size_mean <- 0
  } else {
    clusters_summary <- present_isolates[!grepl("singleton|1pt", names(present_isolates))]
    clusters <- length(clusters_summary)
    cluster_isolates <- sum(clusters_summary)
    cluster_size_median <- stats::median(clusters_summary %>% unlist) %>% round(., 2)
    cluster_size_range <- range(clusters_summary) %>% paste0(collapse = "-")
    cluster_size_mean <- mean(clusters_summary) %>% round(., 2)
  }

  # Characerize absent
  absent <- subset(sum_clustering, grepl("revertant|no feature", names(sum_clustering)))
  no_feature <- sum(absent)
  if (sum(grepl("revertant", names(absent))) >0) {
    revertant_summary <- subset(absent, grepl("revertant", names(absent))) %>% unlist %>% subset(is.na(.) == FALSE)
    revertant_isolates <- sum(revertant_summary)
    revertant_clusters_summary <- subset(revertant_summary, names(revertant_summary) != "revertant_tip")
    revertant_clusters <- length(revertant_clusters_summary)
    if (revertant_clusters == 0) {
      revertant_cluster_size_median <- 0
      revertant_cluster_size_mean <- 0
      revertant_cluster_size_range <- NA
    } else {
      revertant_cluster_size_median <- stats::median(revertant_clusters_summary %>% unlist) %>% round(., 2)
      revertant_cluster_size_mean <- mean(revertant_clusters_summary) %>% round(., 2)
      revertant_cluster_size_range <- range(revertant_clusters_summary) %>% paste0(collapse = "-")
    }
  } else {
    revertant_isolates <- 0
    revertant_clusters <- 0
    revertant_cluster_size_median <- 0
    revertant_cluster_size_mean <- 0
    revertant_cluster_size_range <- NA
  }

  # Phylo frequency
  phylogenetic_events <- sum(singletons + clusters)
  feature_frequency <- (present / num_isolates * 100) %>%  round(., 2)
  phylogenetic_frequency <- (phylogenetic_events / sum(phylogenetic_events + absent) * 100)  %>% round(., 2)
  clustering_frequency <- (clusters / phylogenetic_events * 100) %>% round(., 2)

  results <- cbind.data.frame(present,
                              singletons, singleton_isolates,
                              clusters, cluster_isolates, cluster_size_median, cluster_size_mean, cluster_size_range,
                              no_feature, revertant_isolates, revertant_clusters, revertant_cluster_size_median, revertant_cluster_size_mean, revertant_cluster_size_range,
                              phylogenetic_events, feature_frequency, phylogenetic_frequency, clustering_frequency)
  return(results)
}
