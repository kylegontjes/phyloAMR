#' Generate Ancestral Reconstruction and Clustering Data
#'
#' This function performs ancestral state reconstruction on a given phylogenetic tree and
#' phenotypic data, generates parent-child relationships, and identifies clusters.
#' It returns a list containing the ancestral reconstruction model, parent-child
#' relationships, and tip data merged with clustering information.
#'
#' @param df A dataframe containing isolate numbers and patient IDs columns labelled as `$isolate_no` and `$Patient_ID`, respectively.
#' @param tr A phylogenetic tree object.
#' @param pheno A string specifying the column name in `df` representing the binary phenotype to be analyzed.
#' @param type A character string specifying the type of ancestral reconstruction model (e.g., "discrete").
#' @param model A character string specifying the model used for ancestral state reconstruction (e.g., "ER" for equal rates).
#' @param conf_threshold A numeric value specifying the confidence threshold for filtering parent-child relationships.
#' @param remove_faux A character string ("yes" or "no") indicating whether to remove or relabel faux clusters.
#' @param confidence A character string ("high or "low") specifying the confidence level for identifying the transitions in the clustering data.
#' @param remove_revertant A character string ("yes" or "no") indicating whether to remove revertant calls from the cleaned text string
#' @param collapse_cluster A character string ("yes or "no") indicating whether to create an additional variable where all singletons and cluster isolates (i.e., cluster 1, cluster 2) are collapsed into one a general category (i.e., singleton OR cluster).
#' @return A list containing the following elements:
#'   \describe{
#'     \item{asr_model}{A dataframe of the ancestral state reconstruction results.}
#'     \item{asr_model_statistics}{A dataframe with statistics on the ancestral state reconstruction model}
#'     \item{parent_child_df}{A dataframe of parent-child relationships with additional descriptive information.}
#'     \item{tip_data_df}{A dataframe of tip data merged with clustering information.}
#'   }
#' @export
#' @examples
#' \dontrun{
#' df <- phylosuite::phylosuiteR_data
#' tr <- phylosuite::tr
#' result <- asr_clustering_algorithm(df = df, tr = tr, pheno = "MVB_non_s",
#'     type = "discrete", model = "ER", conf_threshold = 0.875,
#'     remove_faux = "yes", confidence = "high")
#' }

asr_clustering_algorithm <- function(df,tr,pheno,type="discrete",model="ER",conf_threshold=0.875,remove_faux="yes",confidence="high",remove_revertant="yes",collapse_cluster="yes"){
  outcome_str <- df[,pheno] %>% `names<-`(df$isolate_no)
  tr$node.label <- NULL
  # Ancestral character estimation reconstruction modeling of binary outcome
  asr_model <- ape::ace(x=outcome_str,phy=tr,type=type,model=model)
  asr_model$fitER <- cbind(asr_model$lik.anc,node = rownames(asr_model$lik.anc) %>% as.numeric) %>% as.data.frame()
  asr_model_statistics <- get_asr_model_statstics(asr_model)
  # Get Parent Child data
  parent_child_df <- get_parent_child_data(tr=tr,anc_data=asr_model$fitER,pheno_data=outcome_str,conf_threshold = conf_threshold)
  # Annotate parent child data
  parent_child_df <- get_phenotypic_continuation_data(parent_child_df)
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
  results <- list(asr_model=asr_model,
                  asr_model_statistics=asr_model_statistics,
                  parent_child_df = parent_child_df,
                  tip_data_df = tip_data_df)

  return(results)
}

get_asr_model_statstics <- function(asr_model){
  rate1 <- asr_model$rates[1]
  se1<- asr_model$se[1]
  loglik <-  asr_model$loglik
  rate2 <-  asr_model$rates[2]
  se2 <-  asr_model$se[2]
  results <- data.frame(loglik,rate1,se1,rate2,se2)
  return(results)
}

get_parent_child_data <- function(tr,anc_data,pheno_data,conf_threshold=0.875){
  # Tree edge info
  edge <- tr$edge %>% as.data.frame
  colnames(edge) <- c("parent", "child")

  # Prediction data
  anc_data$pred <- ifelse(anc_data$`1` > conf_threshold,1,ifelse(anc_data$`0` >conf_threshold,0,0.5))
  nodes <- anc_data$node %>% unique
  num_edges <- nrow(edge)

  # Dataframe construction
  edge$parent_mle <- rep(NA,num_edges)
  edge$parent_val <- rep(NA,num_edges)
  edge$child_mle <- rep(NA,num_edges)
  edge$child_val <- rep(NA,num_edges)

  # Prediction values for the parent & child
  for(i in nodes){
    edge[edge$parent == i,"parent_mle"] <- anc_data[anc_data$node == i,'1'] *100
    edge[edge$parent == i,"parent_val"] <- anc_data[anc_data$node == i,"pred"]
    edge[edge$child == i,"child_mle"] <- anc_data[anc_data$node == i,'1'] *100
    edge[edge$child == i,"child_val"] <- anc_data[anc_data$node == i,"pred"]
  }

  # Values and name for the tip
  for(i in 1:(min(nodes)-1)){
    edge[edge$child == i,"child_val"] <- pheno_data[[i]]
    edge[edge$child == i,"child_name"] <- names(pheno_data[i])
  }

  return(edge)
}

get_phenotypic_continuation_data <- function(parent_child_df){
  parent_child_df <- parent_child_df %>% dplyr::mutate(transition_any = ifelse(parent_val != child_val,1,0),
                                                transition_high = ifelse(parent_val ==0 & child_val==1 | parent_val==1 & child_val ==0,1,0),
                                                transition_low = ifelse(parent_val ==0.5 & child_val==1 | parent_val==0.5 & child_val ==0,1,0),
                                                gain_low = ifelse(child_val ==1 & parent_val==0.5,1,0),
                                                gain_high  =  ifelse(child_val ==1 & parent_val==0,1,0),
                                                gain_any = ifelse(gain_low==1 | gain_high==1,1,0),
                                                loss_low  =  ifelse(child_val ==0 & parent_val==0.5,1,0) ,
                                                loss_high  =  ifelse(child_val ==0 & parent_val==1,1,0) ,
                                                loss_any = ifelse(loss_low==1 | loss_high==1,1,0),
                                                continuation_any  =  ifelse(parent_val == child_val,1,0),
                                                continuation_high  =  ifelse(parent_val == 1 & child_val ==1 | parent_val == 0 & child_val ==0,1,0),
                                                continuation_low =  ifelse(parent_val == 0.5 & child_val ==0.5,1,0))
  return(parent_child_df)
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
