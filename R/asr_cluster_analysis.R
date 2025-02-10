#' Cluster analysis
#'
#' Description of what the function does.
#'
#' @param x Description of parameter `x`
#' @param y Description of parameter `y`
#' @return Description of return value
#' @export
asr_cluster_analysis <- function(asr_cluster_obj,remove_faux="yes"){
  if(remove_faux =="yes"){
    clustering <- asr_cluster_obj[["tip_data_df"]][["asr_cluster"]] %>% {ifelse(grepl("1pt",.)==T,"r_singleton",.)}
  } else {
    clustering <- asr_cluster_obj[["tip_data_df"]][["asr_cluster"]]
  }

  sum_clustering <-  table(clustering)
  r_isolates <- subset(sum_clustering,grepl("r_",names(sum_clustering)))
  r_ct <- sum(r_isolates)
  r_singleton <- subset(r_isolates,names(r_isolates)=="r_singleton")  %>% ifelse(length(.)==0,0,.)
  r_clusters <-  if(r_singleton==r_ct){0} else { subset(r_isolates,names(r_isolates)!="r_singleton")  %>% unlist %>% subset(is.na(.)==F)}

  r_clusters_num <- sum(r_clusters)
  r_cluster_ct <-  if(r_singleton==r_ct){0} else { subset(r_isolates,names(r_isolates)!="r_singleton")  %>% unlist %>% subset(is.na(.)==F) %>% length}

  r_cluster_median_num <- stats::median(r_clusters %>% unlist) %>% round(.,2)
  r_cluster_range <- range(r_clusters) %>% paste0(collapse = "-")
  r_cluster_mean_num <- mean(r_clusters) %>% round(.,2)

  s <- subset(sum_clustering,!grepl("r_",names(sum_clustering)))
  s_isolates <- sum(s)
  s_rev <- subset(s,grepl("s_rev",names(s))) %>% unlist %>% subset(is.na(.)==F)
  s_rev_ct <- length(s_rev)
  s_rev_num <- sum(s_rev)
  s_rev_median_num <- ifelse(s_rev_num==0,0,stats::median(s_rev %>% unlist))
  s_rev_range <- ifelse(s_rev_num==0,0,range(s_rev) %>% paste0(collapse = "-"))
  s_rev_mean_num <- ifelse(s_rev_num==0,0,mean(s_rev))

  # Phylo frequency
  phylo_events <- sum(r_singleton + r_cluster_ct)
  feature_frequency <- {r_ct / nrow(asr_cluster_obj[["tip_data_df"]]) * 100} %>%  round(.,2)
  phylo_frequency <- {phylo_events / sum(phylo_events + s_isolates) * 100}  %>% round(.,2)
  fixation_frequency <- {r_cluster_ct / phylo_events * 100} %>% round(.,2)

  # Tip data
  num_transitions_any <- asr_cluster_obj[["parent_child_df"]][["transition_any"]] %>% sum
  num_transitions_high <- asr_cluster_obj[["parent_child_df"]][["transition_high"]] %>% sum
  num_loss_any <- asr_cluster_obj[["parent_child_df"]][["loss_any"]] %>% sum
  num_loss_high <- asr_cluster_obj[["parent_child_df"]][["loss_high"]] %>% sum
  num_loss_tip_any <- asr_cluster_obj[["tip_data_df"]][["loss_any"]]%>% sum
  num_loss_tip_high <- asr_cluster_obj[["tip_data_df"]][["loss_high"]]%>% sum
  num_gain_any <- asr_cluster_obj[["parent_child_df"]][["gain_any"]] %>% sum
  num_gain_high <- asr_cluster_obj[["parent_child_df"]][["gain_high"]] %>% sum
  num_gain_tip_any <- asr_cluster_obj[["tip_data_df"]][["gain_any"]]%>% sum
  num_gain_tip_high <- asr_cluster_obj[["tip_data_df"]][["gain_high"]]%>% sum

  results <- cbind.data.frame(r_ct,
                              r_singleton,
                              r_clusters_num,r_cluster_ct,r_cluster_median_num,r_cluster_mean_num,r_cluster_range,
                              s_isolates,s_rev_ct,s_rev_num,s_rev_median_num,s_rev_mean_num,s_rev_range,
                              phylo_events,feature_frequency,phylo_frequency,fixation_frequency,
                              num_transitions_any,
                              num_transitions_high,
                              num_loss_any,num_loss_high,num_loss_tip_any,num_loss_tip_high,
                              num_gain_any,num_gain_high,num_gain_tip_any,num_gain_tip_high)
  return(results)
}
