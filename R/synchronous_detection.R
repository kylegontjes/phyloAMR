#' Synrchronous detection
#'
#' Description of what the function does.
#'
#' @param x Description of parameter `x`
#' @param y Description of parameter `y`
#' @return Description of return value
#' @export
synchronous_detection <- function(comparitor_parent_child_df,trait_parent_child_df){
  gains_t = trait_parent_child_df %>% subset(gain_any==1) %>% .$child
  gains_c = comparitor_parent_child_df  %>% subset(gain_any==1) %>% .$child
  loss_t = trait_parent_child_df %>% subset(loss_any==1) %>% .$child
  loss_c = comparitor_parent_child_df  %>% subset(loss_any==1) %>% .$child

  num_trait_gains =  length(gains_t)
  num_trait_loss = length(loss_t)

  synchronous_gains = gains_t[which(gains_t %in%  gains_c)]
  synchronous_gains_str = paste0(synchronous_gains,collapse = ",")
  synchronous_gains_num = length(synchronous_gains)
  synchronous_gains_prop =   synchronous_gains_num / num_trait_gains
  synchronous_gain_loss = gains_t[which(gains_t %in%  loss_c)]
  synchronous_gain_loss_str = paste0(synchronous_gain_loss,collapse = ",")
  synchronous_gain_loss_num =  length(synchronous_gain_loss)
  synchronous_gain_loss_prop =  synchronous_gain_loss_num / num_trait_gains
  synchronous_loss = loss_t[which(loss_t %in%  loss_c)]
  synchronous_loss_str = paste0(synchronous_loss,collapse = ",")
  synchronous_loss_num = length(synchronous_loss)
  synchronous_loss_prop =  synchronous_loss_num / num_trait_loss
  synchronous_loss_gain = loss_t[which(loss_t %in%  gains_c)]
  synchronous_loss_gain_str = paste0(synchronous_loss_gain,collapse = ",")
  synchronous_loss_gain_num =  length(synchronous_loss_gain)
  synchronous_loss_gain_prop =  synchronous_loss_gain_num / num_trait_loss

  summary <- data.frame(num_trait_gains=num_trait_gains,synchronous_gains=synchronous_gains_str,synchronous_gains_num = synchronous_gains_num,synchronous_gains_prop=synchronous_gains_prop,synchronous_gain_loss=synchronous_gain_loss_str,synchronous_gain_loss_num=synchronous_gain_loss_num,synchronous_gain_loss_prop=synchronous_gain_loss_prop,num_trait_loss=num_trait_loss,synchronous_loss=synchronous_loss_str,synchronous_loss_num=synchronous_loss_num,synchronous_loss_prop=synchronous_loss_prop,synchronous_loss_gain=synchronous_loss_gain_str,synchronous_loss_gain_num=synchronous_loss_gain_num,synchronous_loss_gain_prop=synchronous_loss_gain_prop)

  return(summary)
}
