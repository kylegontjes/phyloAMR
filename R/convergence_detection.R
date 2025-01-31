convergence_detection <- function(comparitor_parent_child_df,trait_parent_child_df){
  gains_t = trait_parent_child_df %>% subset(gain_any==1) %>% .$child
  gains_c = comparitor_parent_child_df  %>% subset(gain_any==1) %>% .$child
  loss_t = trait_parent_child_df %>% subset(loss_any==1) %>% .$child
  loss_c = comparitor_parent_child_df  %>% subset(loss_any==1) %>% .$child

  num_trait_gains =  length(gains_t)
  num_trait_loss = length(loss_t)

  convergent_gains = gains_t[which(gains_t %in%  gains_c)]
  convergent_gains_num = length(convergent_gains)
  convergent_gains_prop =   convergent_gains_num / num_trait_gains
  convergent_loss = loss_t[which(loss_t %in%  loss_c)]
  convergent_loss_num = length(convergent_loss)
  convergent_loss_prop =  convergent_loss_num / num_trait_loss

  summary <- data.frame(num_trait_gains=num_trait_gains,convergent_gains_num = convergent_gains_num,convergent_gains_prop=convergent_gains_prop,num_trait_loss=num_trait_loss,convergent_loss_num=convergent_loss_num,convergent_loss_prop=convergent_loss_prop)

  results = list(summary=summary,
                 convergent_gains=convergent_gains,
                 convergent_loss=convergent_loss
  )

  return(results)
}
