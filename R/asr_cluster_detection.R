#' Cluster detection
#'
#' Description of what the function does.
#'
#' @param df Dataframe with tip name variable and phenotype
#' @param tr Phylogenetic tree
#' @param tip_name_var Name of variable containing tip names in df
#' @param patient_id Name of variable containing patient IDs, can be combined with faux_clusters option to factor into whether a cluster should have >1 patient. (Optional)
#' @param pheno Name of phenotype variable in df
#' @param parent_child_df Parent child dataframe from asr() object
#' @param node_states Whether the reconstruction was "joint" or "marginal"
#' @param confidence Whether to use 'high' (i.e., 0 -> 1) or 'low' (0 -> 0.5) confidence transitions when determining clustering. ONLY USAGE FOR MARGINAL STATE RECONSTRUCTIONS
#' @param faux_clusters Two options exist: relabel (i.e., specify clusters as having 1pt) or remove (i.e., consider these singletons) (Optional)
#' @param remove_revertant Whether to remove revertant episodes from consideration in a cleaned
#' @param collapse_cluster Whether to create a variable that collapses cluster calls into one category
#' @return A tip-only dataframe with inferences on the history of these strains. Can be merged with parent_child_df from asr() if desired
#' @export
asr_cluster_detection <- function(df,tr,tip_name_var,patient_id=NULL,pheno,parent_child_df,node_states="joint",confidence=NULL,faux_clusters=F,remove_revertant="yes",collapse_cluster="yes"){
  # Check if states are as desired
  check_joint_confidence(node_states,confidence)
  # Check faux_cluster
  check_faux_clusters(patient_id,faux_clusters)

  # Set up data frame
  clustering_data <- data.frame(isolate_no = df[[tip_name_var]])
  rownames(clustering_data) <- clustering_data$isolate_no
  # If patient ID is not null
  if(is.null(patient_id)==F){
    clustering_data$patient_id = df[[patient_id]]
  }

  # Get clustering data
  clustering_data$asr_cluster <- sapply(clustering_data$isolate_no,FUN=get_clustering_data,parent_child_df,node_states,confidence)

  # Clean clusters
  if(faux_clusters == "remove"){
    clustering_data$asr_cluster <- remove_faux_clusters(clustering_data)
  }

  if(faux_clusters == "relabel"){
    clustering_data$asr_cluster <- remove_faux_clusters(clustering_data)
  }

  # Simplify strings
  clustering_data$asr_cluster_renamed <- simplify_acr_clustering_string(string = "asr_cluster",df = clustering_data,tr = tr,faux_clusters =faux_clusters,remove_revertant = remove_revertant)

  # Simplify to cluster
  if(collapse_cluster=="yes"){
    clustering_data$asr_cluster_collapsed <- group_by_category(clustering_data$asr_cluster)
  }

  # Get tip data and merge it with the clustering data
  tip_data <- parent_child_df %>% subset(child <= ape::Nnode(tr)+1) %>% dplyr::mutate(isolate_no=child_name)
  tip_data_df <- tip_data %>% dplyr::left_join(.,clustering_data) %>% `rownames<-`(.$isolate_no) %>% select(-isolate_no)

  # We don't need the parent_child dataframe here, so we are only providing the tip based data
  results <-tip_data_df
  return(results)
}

check_joint_confidence <- function(node_states,confidence){
  if(node_states == "joint" & is.null(confidence)== F
     ){
    stop("When using joint reconstruction states, you must specify confidence as NULL")
  }
}

check_faux_clusters <- function(patient_id,faux_clusters){
  if(is.null(patient_id) == T & faux_clusters != F){
    stop("Must have patient id variable for when faux clusters is specified")
  }
  if(!faux_clusters %in% c(F,"relabel","remove")){
    stop("Must specify faux cluster under options of F, relabel, or remove")
  }
}


get_clustering_data <- function (isolate, edge_data,node_states, confidence){
  tip_row <- edge_data %>% subset(child_name == isolate)
  classification <- c()
  if (tip_row$child_val == 0) {
    if(node_states == "joint"){
      classification <- ifelse(tip_row[,"loss"] == 1, "revertant_tip",
                               "no feature")
    } else {
      classification <- ifelse(tip_row[,paste0("loss_",confidence)] == 1, "revertant_tip",
                               "no feature")
    }
  }

  if (tip_row$child_val == 1) {
    if(node_states == "joint"){
      classification <- ifelse(tip_row[,"gain"] == 1, "singleton",
                               paste0("cluster_", get_largest_anc_cluster(tip_row, edge_data) %>% {ifelse(.=="root","root",.)}))
    } else {

      classification <- ifelse(tip_row[,paste0("gain_",confidence)] == 1, "singleton",
                               paste0("cluster_", get_largest_anc_cluster(tip_row, edge_data) %>% {ifelse(.=="root","root",.)}))

    }
  }
  if (classification == "singleton") {
    classification <- reclassify_singletons(node_data = tip_row,
                                                         edge_data = edge_data,node_status=node_states,confidence=confidence)
  }
  if (classification == "no feature") {
    classification <- get_largest_anc_reversion(node_data = tip_row,
                                                                 edge_data = edge_data,node_states=node_states,confidence=confidence)

    classification <- ifelse(classification == "no feature",classification,paste0("revertant_",classification))
  }
  return(classification)
}

get_largest_anc_cluster <- function(node_data,edge_data){
  output <- node_data
  repeat{
    output <- get_parent_node(output, edge_data)
    if(nrow(output)==0){
      output <- c("root")
      break
    }
    else {
      if(node_states =="joint"){
        if(output[,"parent_val"]==0){
          output <- output$child
          break
        }
      }
      if(node_states == "marginal"){
        if(confidence=="high"){
          if(output[,"parent_val"]==0){
            output <- output$child
            break
          }
        } else {
          if(output[,"parent_val"]==0.5){
            output <- output$child
            break
          }
        }
      }
    }
  }
  return(output)
}

reclassify_singletons <- function (node_data, edge_data,node_status,confidence){
  previous <- get_parent_node(node_data, edge_data)
  if(node_status == "joint"){
    classification <- ifelse(previous[,"loss"] == 0, "singleton",
                             get_largest_anc_cluster(previous, edge_data))
  } else {
    classification <- ifelse(previous[,paste0("loss_",confidence)] == 0, "singleton",
                             get_largest_anc_cluster(previous, edge_data))
  }
  return(classification)
}

get_largest_anc_reversion <- function (node_data, edge_data,node_states,confidence){
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
    else
      if(node_states =="joint"){
        if(output[,"loss"]==1){
          output <- output$child
          break
        }
      }
      if(node_states == "marginal"){
        if (output[,paste0("loss_",confidence)] == 1) {
          output <- output$child
          break
        }
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
  cluster_size <- table(cluster_data$asr_cluster) %>% subset(grepl("r_",names(.))) %>% subset(!grepl("singleton",names(.)))
  # Isolates where gain at internal node was called, but no other isolate was connected to it
  true_clusters <- subset(cluster_size,cluster_size>1) %>% names
  not_cluster <- names(cluster_size[!names(cluster_size) %in% true_clusters])
  # Clusters with only one patient
  cluster_size_pts <- sapply(true_clusters,FUN=function(x){
    cluster_data  %>% subset(asr_cluster == x) %>% .$patient_id %>% unique %>% length
  })
  not_more_than_one_pt <- subset(cluster_size_pts,cluster_size_pts==1) %>% names(.)
  # Reclassify these isolates as singletons
  cluster_fin <-  ifelse(cluster_data$asr_cluster  %in% c(not_more_than_one_pt,not_cluster),"singleton",cluster_data$asr_cluster)

  return(cluster_fin)
}

relabel_faux_clusters <- function(cluster_data) {
  cluster_size <- table(cluster_data$asr_cluster) %>% subset(grepl("r_",
                                                                   names(.))) %>% subset(!grepl("singleton", names(.)))

  # Isolates where gain at internal node was called, but no other isolate was connected to it
  true_clusters <- subset(cluster_size, cluster_size > 1) %>%
    names
  not_cluster <- names(cluster_size[!names(cluster_size) %in%
                                      true_clusters])

  # Clusters with only one patient
  cluster_size_pts <- sapply(true_clusters, FUN = function(x) {
    cluster_data %>% subset(asr_cluster == x) %>% .$patient_id %>%
      unique %>% length
  })
  not_more_than_one_pt <- subset(cluster_size_pts, cluster_size_pts ==
                                   1) %>% names(.)
  # Rename only one isolate as singleton & add _1pt_only label to the faux single patient cluster
  cluster_fin <- ifelse(cluster_data$asr_cluster %in% not_cluster,"singleton",
                        ifelse(cluster_data$asr_cluster %in% not_more_than_one_pt,paste0(cluster_data$asr_cluster,"_1pt_only"),
                               cluster_data$asr_cluster))

  return(cluster_fin)
}

simplify_acr_clustering_string <-  function(string,df,tr,faux_clusters,remove_revertant){
  # add rownames
  rownames(df) <- df$isolate_no
  #Reorder vector to tree plotting
  df_ordered <- df %>% dplyr::select(isolate_no,paste0(string)) %>% .[match(ggtree::get_taxa_name(ggtree::ggtree(tr)),df$isolate_no),]

  # Convert Names to Easy to Report String
  df_ordered$string <- df_ordered %>% dplyr::select(paste0(string)) %>% unlist %>%  {ifelse(. == "singleton","Singleton",
                                                                                            ifelse(. == "no feature","No Feature",
                                                                                                   ifelse(.!="revertant_tip" & grepl("revertant",.),"No Feature",
                                                                                                          ifelse(.=="revertant_tip","Revertant Tip",.))))}
  # Name clusters
  num_clusters <- length(df_ordered$string  %>% subset(.,!. %in% c("No Feature","Singleton","Reversion Tip")) %>% subset(.,grepl("1pt",.)==F) %>% unique)
  if(num_clusters >0){
    cluster_names <- df_ordered$string  %>% subset(.,!. %in% c("No Feature","Singleton","Reversion Tip")) %>% subset(.,grepl("1pt",.)==F) %>% unique %>% `names<-`(paste0("Cluster ",1:length(.))) %>% stats::setNames(names(.), .)
    df_ordered$string <- dplyr::recode(df_ordered$string, !!! cluster_names)

    if(faux_clusters=="rename"){
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
         ifelse(string == "singleton","singleton",
                ifelse(grepl("cluster_",string)==T,"cluster",
                       ifelse(grepl("revertant_tip",string)==T,"revertant tip",
                              ifelse(grepl("revertant_",string)==T,"revertant lineage",NA)))))
}
