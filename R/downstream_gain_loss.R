#' Downstream gain loss
#'
#' Description of what the function does.
#'
#' @param x Description of parameter `x`
#' @param y Description of parameter `y`
#' @return Description of return value
#' @export
downstream_gain_loss <- function(comparitor_parent_child_df,parent_child_df,tr){
  stretches = get_trait_traces_on_tree(parent_child_df,tr)
  downstream_changes = get_gain_loss_on_stretches(comparitor_parent_child_df = comparitor_parent_child_df,stretches)
  return(downstream_changes)
}

get_trait_traces_on_tree <- function(parent_child_df,tr){
  all_possible_paths = nodepath(tr)

  gains = parent_child_df %>% subset(gain_any==1) %>% .$child
  tip_gains = subset(gains,gains < min(tr$edge[1,]))
  ancestral_gains = subset(gains,!gains %in% tip_gains)

  gain_paths = lapply(all_possible_paths,FUN=function(x,gains){
    didit = x %in% gains
    if(sum(didit)==0){
      string=NA
    } else {
      didit_pos = which(didit)
      string=x[didit_pos:length(x)]
    }

    return( string)
  },gains=ancestral_gains) %>% subset(is.na(.)==F)

  trait_found = unique(c(subset(parent_child_df,child_val==1) %>% .$child,subset(parent_child_df,parent_val==1) %>% .$parent))

  paths_w_trait = lapply(gain_paths,FUN=function(x,trait_found){
    results = x %in% trait_found
    names(results) = x
    return(results)
  },trait_found=trait_found)

  find_true_stretches <- function(x) {
    rle_result <- rle(x)  # Run-length encoding
    nodes = names(x)
    starts <- cumsum(c(1, head(rle_result$lengths, -1)))  # Start positions
    stretches <- data.frame(
      startnode =names( starts[rle_result$values]),start = starts[rle_result$values],
                            end = starts[rle_result$values] + rle_result$lengths[rle_result$values] - 1)
    paths = list()
    if(nrow(stretches)>0){

      for(i in 1:nrow(stretches)){
        start= stretches[i,"start"]
        end = stretches[i,"end"]
        paths[[i]] <- as.vector(nodes[start:end])
      }
    }
    paths = unlist(paths)


    return(paths)
  }

  true_stretches = lapply(paths_w_trait,find_true_stretches)
  return(true_stretches)
}


get_gain_loss_on_stretches = function(comparitor_parent_child_df, stretches){
  gains_c = comparitor_parent_child_df  %>% subset(gain_any==1) %>% .$child
  loss_c = comparitor_parent_child_df  %>% subset(loss_any==1) %>% .$child

  get_downstream_nodes = function(stretches){
    lapply(stretches,FUN=function(x){x[-1]})
  }

  downstream = get_downstream_nodes(stretches)
  downstream_positions = stretches %>% unlist %>% unique %>% sort

  num_downstream = length(downstream_positions)
  gains = downstream_positions[which(downstream_positions %in%  gains_c)]
  gains_str = paste0(gains,collapse=",")
  gains_num = length(gains)
  gains_prop = gains_num / num_downstream
  loss =  downstream_positions[which(downstream_positions %in%  loss_c)]
  loss_str = paste0(loss,collapse=",")
  loss_num =length(loss)
  loss_prop = loss_num / num_downstream

  # Collated lists
  get_gain_paths = function(stretches){
    gains =  lapply(stretches,FUN=function(x){x[1]}) %>%unlist %>% unique %>% sort
    # group by stretches with same start
    paths_from_gain = lapply(gains,FUN=function(gain){
      stretches[lapply(stretches,FUN=function(y){gain %in% y})  %>% unlist %>% which] %>% unlist %>% as.numeric  %>% sort(decreasing = T) %>% unique()
    })

    downstream_merged_paths = lapply(paths_from_gain,FUN=function(x){
      without_downstream_but_merged =          subset(x,!x%in%gains)

      return(without_downstream_but_merged)})
    names(downstream_merged_paths) = lapply(paths_from_gain,FUN=function(x){
      subset(x,x%in%gains)
    }) %>% unlist

    return(downstream_merged_paths)
  }
  merged_paths = get_gain_paths(stretches = stretches )
  # Num paths
  num_stretches = length(merged_paths)
  #Gains
  stretches_w_gains = lapply(merged_paths,FUN=function(x){
    x[which(x %in% gains_c)] %>% sum %>% {ifelse(.>0,T,F)}
  }) %>% unlist %>% which %>% names %>% as.numeric %>% {ifelse(length(.>0),.,NA)}

  stretches_w_gains_num = length(stretches_w_gains)
  stretches_w_gains_prop  = stretches_w_gains_num / num_stretches
  # Losses
  stretches_w_losses = lapply(merged_paths,FUN=function(x){
    x[which(x %in% loss_c)] %>% sum %>% {ifelse(.>0,T,F)}
  })  %>% unlist %>% which %>% names %>% as.numeric %>% {ifelse(length(.>0),.,NA)}
  stretches_w_losses_num = length(stretches_w_losses)
  stretches_w_losses_prop = stretches_w_losses_num / num_stretches

  summary = data.frame(num_stretches=num_stretches,stretches_w_gains = stretches_w_gains,stretches_w_gains_num=stretches_w_gains_num,stretches_w_gains_prop=stretches_w_gains_prop,stretches_w_losses = stretches_w_losses,stretches_w_losses_num,stretches_w_losses_prop=stretches_w_losses_prop,num_downstream=num_downstream,gains=gains_str,gains_num=gains_num,gains_prop=gains_prop,loss=loss_str,loss_num=loss_num,loss_prop=loss_prop)
  return(summary)
}

