#' asr_cluster_analysis: Cluster analysis
#'
#' This function calculates the clustering statistics
#'
#' @param tip_data_df Tip data frame object from asr_cluster_detection()
#' @param remove_faux whether to remove faux isolates or not
#' @return A dataframe with data on the phylogenetics of a trait
#' @export
asr_cluster_analysis <- function(tip_data_df, remove_faux = "yes") {
  clustering <- tip_data_df[["asr_cluster"]]

  if (remove_faux == "yes") {
    clustering <- ifelse(grepl("1pt", clustering), "singleton", clustering)
  }

  sum_clustering <-  table(clustering)
  num_isolates  <-  sum(sum_clustering)
  # Characterize present
  present_isolates <- subset(sum_clustering, grepl("cluster|singleton", names(sum_clustering)))
  present <- sum(present_isolates)
  singletons <- subset(present_isolates, names(present_isolates) == "singleton")
  singletons <- ifelse(length(singletons) == 0, 0, singletons)
  if (singletons == present) {
    clusters <- 0
    cluster_isolates <-  0
    cluster_size_median <- 0
    cluster_size_range <- NA
    cluster_size_mean <- 0
  } else {
    clusters_summary <-  subset(present_isolates, names(present_isolates) != "singleton")  %>% unlist %>% subset(is.na(.) == FALSE)
    clusters <- length(clusters_summary)
    cluster_isolates <- sum(clusters_summary)
    cluster_size_median <- stats::median(clusters_summary %>% unlist) %>% round(., 2)
    cluster_size_range <- range(clusters_summary) %>% paste0(collapse = "-")
    cluster_size_mean <- mean(clusters_summary) %>% round(., 2)
  }

  #Characerize absent
  absent <- subset(sum_clustering, !grepl("cluster|singleton", names(sum_clustering)))
  no_feature <- sum(absent)
  if (sum(grepl("revertant", names(absent))) >0) {
    revertant_summary <- subset(absent, grepl("revertant", names(absent))) %>% unlist %>% subset(is.na(.) == FALSE)
    revertant_isolates <- sum(revertant_summary)
    revertant_lineages_summary <- subset(revertant_summary, names(revertant_summary) != "revertant_tip")
    revertant_lineages <- length(revertant_lineages_summary)
    if (revertant_lineages == 0) {
      revertant_lineage_size_median <- 0
      revertant_lineage_size_mean <- 0
      revertant_lineage_size_range <- NA
    } else {
      revertant_lineage_size_median <- stats::median(revertant_lineages_summary %>% unlist) %>% round(., 2)
      revertant_lineage_size_mean <- mean(revertant_lineages_summary) %>% round(., 2)
      revertant_lineage_size_range <- range(revertant_lineages_summary) %>% paste0(collapse = "-")
    }
  } else {
    revertant_isolates <- 0
    revertant_lineages <- 0
    revertant_lineage_size_median <- 0
    revertant_lineage_size_mean <- 0
    revertant_lineage_size_range <- NA
  }

  # Phylo frequency
  phylogenetic_events <- sum(singletons + clusters)
  feature_frequency <- (present / num_isolates * 100) %>%  round(., 2)
  phylogenetic_frequency <- (phylogenetic_events / sum(phylogenetic_events + absent) * 100)  %>% round(., 2)
  clustering_frequency <- (clusters / phylogenetic_events * 100) %>% round(., 2)

  results <- cbind.data.frame(present,
                              singletons,
                              clusters, cluster_isolates, cluster_size_median, cluster_size_mean, cluster_size_range,
                              no_feature, revertant_isolates, revertant_lineages, revertant_lineage_size_median, revertant_lineage_size_mean, revertant_lineage_size_range,
                              phylogenetic_events, feature_frequency, phylogenetic_frequency, clustering_frequency)
  return(results)
}
