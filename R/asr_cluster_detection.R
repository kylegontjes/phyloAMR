#' Cluster detection
#'
#' Description of what the function does.
#'
#' @param x Description of parameter `x`
#' @param y Description of parameter `y`
#' @return Description of return value
#' @export
asr_cluster_detection <- function(df,tr,pheno,parent_child_df,remove_faux="yes",confidence="high",remove_revertant="yes",collapse_cluster="yes"){
  # Get clustering data
  clustering_data <- cbind.data.frame(isolate_no = df$isolate_no,
                                      Patient_ID = df$Patient_ID,
                                      asr_cluster = lapply(df$isolate,FUN=get_clustering_data,parent_child_df,confidence) %>% unlist)
# Clean clusters
if(remove_faux == "yes"){
  clustering_data$asr_cluster <- remove_faux_clusters(clustering_data)
} else {
  clustering_data$asr_cluster <- relabel_faux_clusters(clustering_data)
}

# Simplify strings
clustering_data$asr_cluster_renamed <- simplify_acr_clustering_string(string = "asr_cluster",df = clustering_data,tr = tr,remove_faux =remove_faux,remove_revertant = remove_revertant)
# Simplify to cluster
if(collapse_cluster=="yes"){
  clustering_data$asr_cluster_collapsed <- group_by_category(clustering_data$asr_cluster)
}

# Get tip data and merge it with the clustering data
tip_data <- parent_child_df %>% subset(child <= ape::Nnode(tr)+1) %>% dplyr::mutate(isolate_no=child_name)
tip_data_df <- tip_data %>% dplyr::left_join(.,clustering_data) %>% `rownames<-`(.$isolate_no)

# List containing all of the results from this pipeline
results <- list(parent_child_df = parent_child_df,
                tip_data_df = tip_data_df)
return(results)
}

get_clustering_data <- function (isolate, edge_data, confidence)
{
  tip_row <- edge_data %>% subset(child_name == isolate)
  classification <- c()
  if (tip_row$child_val == 0) {
    classification <- ifelse(tip_row[,paste0("loss_",confidence)] == 1, "s_rev_tip",
                             "no feature")
  }
  if (tip_row$child_val == 1) {
    classification <- ifelse(tip_row[,paste0("gain_",confidence)] == 1, "singleton",
                             paste0("r_", get_largest_anc_cluster(tip_row, edge_data) %>%
                                      .$child))
  }
  if (classification == "singleton") {
    classification <- paste0("r_", reclassify_singletons(node_data = tip_row,
                                                         edge_data = edge_data,confidence=confidence))
  }
  if (classification == "no feature") {
    classification <- paste0("s_rev_", get_largest_anc_reversion(node_data = tip_row,
                                                                 edge_data = edge_data,confidence=confidence)) %>% {
                                                                   ifelse(. == "s_rev_no feature", "no feature", .)
                                                                 }
  }
  return(classification)
}

get_largest_anc_cluster <- function(node_data,edge_data){
  output <- node_data
  repeat{
    output <- get_parent_node(output, edge_data)
    if(nrow(output)==0){
      output <- node_data
      repeat{
        output <- get_parent_node(output, edge_data)
        if(nrow(output)==0){
          output <- c("no_0.5_in_path")
          break
        } else  {if (output$parent_val == 0.5){
          break
        }
        }
      }
      break
    } else  {if (output$parent_val == 0) {
      break
    }
    }
  }
  return(output)
}

reclassify_singletons <- function (node_data, edge_data,confidence){
  previous <- get_parent_node(node_data, edge_data)
  classification <- ifelse(previous[,paste0("loss_",confidence)] == 0, "singleton",
                           get_largest_anc_cluster(previous, edge_data) %>% .$child)
  return(classification)
}



get_largest_anc_reversion <- function (node_data, edge_data,confidence){
  output <- node_data
  repeat {
    output <- get_parent_node(output, edge_data)
    if (nrow(output) == 0) {
      output <- "no feature"
      break
    }
    else if (output$parent_val == 1) {
      output <- output$child
      break
    }
    else if (output[,paste0("loss_",confidence)] == 1) {
      output <- output$child
      break
    }
  }
  return(output)
}

get_parent_node <- function(node_data,edge_data){
  parent_node <- node_data$parent
  step_up <- subset(edge_data,child==parent_node)
  return(step_up)
}


remove_faux_clusters <- function(cluster_data){
  cluster_size <- table(cluster_data$asr_cluster) %>% subset(grepl("r_",names(.))) %>% subset(!grepl("r_singleton",names(.)))
  # Isolates where gain at internal node was called, but no other isolate was connected to it
  true_clusters <- subset(cluster_size,cluster_size>1) %>% names
  not_cluster <- names(cluster_size[!names(cluster_size) %in% true_clusters])
  # Clusters with only one patient
  cluster_size_pts <- sapply(true_clusters,FUN=function(x){
    cluster_data  %>% subset(asr_cluster == x) %>% .$Patient_ID %>% unique %>% length
  })
  not_more_than_one_pt <- subset(cluster_size_pts,cluster_size_pts==1) %>% names(.)
  # Reclassify these isolates as singletons
  cluster_fin <-  ifelse(cluster_data$asr_cluster  %in% c(not_more_than_one_pt,not_cluster),"r_singleton",cluster_data$asr_cluster)

  return(cluster_fin)
}

relabel_faux_clusters <- function(cluster_data) {
  cluster_size <- table(cluster_data$asr_cluster) %>% subset(grepl("r_",
                                                                   names(.))) %>% subset(!grepl("r_singleton", names(.)))

  # Isolates where gain at internal node was called, but no other isolate was connected to it
  true_clusters <- subset(cluster_size, cluster_size > 1) %>%
    names
  not_cluster <- names(cluster_size[!names(cluster_size) %in%
                                      true_clusters])

  # Clusters with only one patient
  cluster_size_pts <- sapply(true_clusters, FUN = function(x) {
    cluster_data %>% subset(asr_cluster == x) %>% .$Patient_ID %>%
      unique %>% length
  })
  not_more_than_one_pt <- subset(cluster_size_pts, cluster_size_pts ==
                                   1) %>% names(.)
  # Rename only one isolate as singleton & add _1pt_only label to the faux single patient cluster
  cluster_fin <- ifelse(cluster_data$asr_cluster %in% not_cluster,"r_singleton",
                        ifelse(cluster_data$asr_cluster %in% not_more_than_one_pt,paste0(cluster_data$asr_cluster,"_1pt_only"),
                               cluster_data$asr_cluster))

  return(cluster_fin)
}

simplify_acr_clustering_string <-  function(string,df,tr,remove_faux,remove_revertant){
  # add rownames
  rownames(df) <- df$isolate_no
  #Reorder vector to tree plotting
  df_ordered <- df %>% dplyr::select(isolate_no,paste0(string)) %>% .[match(ggtree::get_taxa_name(ggtree::ggtree(tr)),df$isolate_no),]

  # Convert Names to Easy to Report String
  df_ordered$string <- df_ordered %>% dplyr::select(paste0(string)) %>% unlist %>%  {ifelse(. == "r_singleton","Singleton",
                                                                                            ifelse(. == "no feature","No Feature",
                                                                                                   ifelse(.!="s_rev_tip" & grepl("s_rev_",.),"No Feature",
                                                                                                          ifelse(.=="s_rev_tip","Reversion Tip",.))))}
  # Name clusters
  num_clusters <- length(df_ordered$string  %>% subset(.,!. %in% c("No Feature","Singleton","Reversion Tip")) %>% subset(.,grepl("1pt",.)==F) %>% unique)
  if(num_clusters >0){
    cluster_names <- df_ordered$string  %>% subset(.,!. %in% c("No Feature","Singleton","Reversion Tip")) %>% subset(.,grepl("1pt",.)==F) %>% unique %>% `names<-`(paste0("Cluster ",1:length(.))) %>% stats::setNames(names(.), .)
    df_ordered$string <- dplyr::recode(df_ordered$string, !!! cluster_names)

    if(remove_faux=="no"){
      # Name faux clusters
      faux_cluster_names <- df_ordered$string %>% subset(.,grepl("1pt",.)==T) %>% unique %>% `names<-`(paste0("Single Patient Cluster ",1:length(.))) %>% stats::setNames(names(.), .)
      df_ordered$string <- dplyr::recode(df_ordered$string, !!! faux_cluster_names)
      df_ordered$string <- factor(levels = c("No Feature","Reversion Tip","Singleton",paste0("Cluster ",1:length(cluster_names)),paste0("Single Patient Cluster ",1:length(faux_cluster_names))),x = df_ordered$string)
    } else {
      # Remove single patient calls
      df_ordered$string <- ifelse(grepl("1pt",df_ordered$string)==T,"Singleton",df_ordered$string)
      # Factor String
      df_ordered$string <- factor(levels = c("No Feature","Reversion Tip","Singleton",paste0("Cluster ",1:length(cluster_names))),x = df_ordered$string)
    }
  }
  if(remove_revertant=="yes"){
    df_ordered$string <- dplyr::recode(df_ordered$string,"Reversion Tip" = "No Feature")
  }

  rownames(df_ordered) <- df_ordered$isolate_no
  # Match order back to original dataframe
  df_reset_order <- df_ordered[match(rownames(df),rownames(df_ordered)),]
  return(df_reset_order$string)
}

group_by_category <- function(string){
  ifelse(string =="no feature","no feature",
         ifelse(string == "r_singleton","r singleton",
                ifelse(grepl("r_",string)==T,"r cluster",
                       ifelse(grepl("s_rev_tip",string)==T,"s revertant tip",
                              ifelse(grepl("s_rev_",string)==T,"s revertant lineage",NA)))))
}
