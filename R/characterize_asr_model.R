#' Characterize ancestral state reconstruction model
#'
#' Function to report model statistics for a binary ancestral state reconstruction model
#'
#' @param corHMM_obj Ancestral reconstruction model from corHMM::corHMM
#' @return Dataframe containing model fit (i.e., log-likelihood and AIC), rates, and model type (i.e., ER or ARD)#'   \describe{
#'     \item{model}{Model}
#'     \item{number_parameters}{Number of parameters in model}
#'     \item{number_rate_categories}{Number of rate categories}
#'     \item{rate1}{Rate 1}
#'     \item{rate2}{Rate 2}
#'     \item{loglik}{Log-likelihood of model}
#'     \item{AIC}{Akaike information criteria}
#'     \item{AICc}{Sample size corected Akaike information criteria}
#'   }
#' @export
characterize_asr_model <- function(corHMM_obj) {
  # Identify the model
  ## Number of parameters in model
  number_parameters <- max(unlist(corHMM_obj$index.mat), na.rm = TRUE)
  ## Identify the model
  model <- ifelse(number_parameters == 1, "ER", "ARD")

  # Number of rate categories
  number_rate_categories <- corHMM_obj$rate.cat

  # Capture the rates
  rate1 <- round(corHMM_obj$solution[1, 2], 2)
  rate2 <- ifelse(model == "ARD", round(corHMM_obj$solution[2, 1], 2), NA)

  # Log likelihood
  loglik <- corHMM_obj$loglik

  # Akaike information criterium
  ## Raw output
  AIC <- corHMM_obj$AIC
  ## Sample-size corrected
  AICc <- corHMM_obj$AICc

  # Data output
  data <- data.frame(model, number_parameters, number_rate_categories, rate1, rate2, loglik, AIC, AICc)

  return(data)
}
