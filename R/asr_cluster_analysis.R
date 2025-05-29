#' Analysis of trait clustering on a phylogeny
#'
#' This function calculates the clustering statistics (e.g., number of singletons, clusters, and no features).
#' Alongside the frequency of these events, the size of resistance and revertant clusters are reported.
#' Two additional statistsics are employed:
#' 1. Phylogenetic frequency: No. clusters + singletons / total opportunities (i.e., phylogenetic events + tips without the trait)
#' 2. Clustering frequency:  No. clusters / total phylogenetic events (i.e., singleton events + clusters)
#'
#' @param tip_data_df Tip data frame object from asr_cluster_detection()
#' @return Dataframe row with data on the phylogenetic clustering of the trait.
#' @export
asr_cluster_analysis <- function(tip_data_df) {
  # Check if it has the expected contents
  check_asr_content(tip_data_df)

  # Clustering string
  clustering <- tip_data_df[["asr_cluster"]]

  # Summary data on number of isolates and membership in clusters
  sum_clustering <-  table(clustering)
  num_isolates  <-  length(clustering)

  # Characterize present
  present_isolates <- sum_clustering[grepl("cluster|singleton", names(sum_clustering)) & !grepl("revertant", names(sum_clustering))]
  present <- sum(present_isolates)
  ## Singleton calls
  singletons <- sum(present_isolates[grepl("singleton", names(present_isolates))], sum(grepl("1pt", names(present_isolates))))
  singleton_isolates <- sum(present_isolates[grepl("singleton|1pt", names(present_isolates))])
  ## Cluster calls
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
    cluster_size_median <- round(stats::median(unlist(clusters_summary)), 2)
    cluster_size_range <- paste0(range(clusters_summary), collapse = "-")
    cluster_size_mean <- round(mean(clusters_summary), 2)
  }

  # Characerize absent
  absent <- subset(sum_clustering, grepl("revertant|no feature", names(sum_clustering)))
  no_feature <- sum(absent)
  ## Revertants
  if (sum(grepl("revertant", names(absent))) > 0) {
    revertant_summary <- subset(absent, grepl("revertant", names(absent)))
    revertant_isolates <- sum(revertant_summary)
    revertant_clusters_summary <- subset(revertant_summary, names(revertant_summary) != "revertant_tip")
    revertant_clusters <- length(revertant_clusters_summary)
    if (revertant_clusters == 0) {
      revertant_cluster_size_median <- 0
      revertant_cluster_size_mean <- 0
      revertant_cluster_size_range <- 0
    } else {
      revertant_cluster_size_median <- round(stats::median(unlist(revertant_clusters_summary)), 2)
      revertant_cluster_size_mean <- round(mean(revertant_clusters_summary), 2)
      revertant_cluster_size_range <- paste0(range(revertant_clusters_summary), collapse = "-")
    }
  } else {
      revertant_isolates <- 0
      revertant_clusters <- 0
      revertant_cluster_size_median <- 0
      revertant_cluster_size_mean <- 0
      revertant_cluster_size_range <- NA
    }

  # Frequency statistics
  feature_frequency <- round(present / num_isolates * 100, 2)
  ## Phylogenetic frequency: No. clusters + singletons / total opportunities (i.e., phylogenetic events + tips without the trait)
  phylogenetic_events <- sum(singletons + clusters)
  phylogenetic_frequency <- round(phylogenetic_events / sum(phylogenetic_events + absent) * 100, 2)
  ## Clustering frequency: No. clusters / total phylogenetic events (i.e., singleton events + clusters)
  clustering_frequency <- round(clusters / phylogenetic_events * 100, 2)

  results <- cbind.data.frame(present,
                              singletons, singleton_isolates,
                              clusters, cluster_isolates, cluster_size_median, cluster_size_mean, cluster_size_range,
                              no_feature, revertant_isolates, revertant_clusters, revertant_cluster_size_median, revertant_cluster_size_mean, revertant_cluster_size_range,
                              phylogenetic_events, feature_frequency, phylogenetic_frequency, clustering_frequency)
  return(results)
}
