# Univariable regression
univariable_regression <- function (variable, outcome, data){
  dataset <- data %>% select(outcome, paste0(variable))
  model <- paste(outcome, paste(variable, collapse = " + "), sep = "~")
  output <- glm(model, data = dataset, family = "binomial")
  ci <- suppressMessages(confint(output))
  final <- cbind(exp(cbind(OR = coef(output), ci)) %>% round(., 2),
                 abs(summary(output)$coefficients[, "Pr(>|z|)"]) %>%
                   round(., 4)) %>% subset(rownames(.) != "(Intercept)") %>%
    `colnames<-`(c("OR", "2.5%", "97.5%", "p_value")) %>%
    as.data.frame %>% mutate(`OR (95% CI)` = paste0(OR, " (", `2.5%`, "-", `97.5%`, ")")) %>% select(`OR (95% CI)`, p_value)
  return(final)
}

univariable_regression_table <- function (outcome, dataset, variables){
  datatable <- lapply(
    variables,
    FUN = function(x) {
      datatable <- univariable_regression(x, outcome, dataset) %>%
        as.data.frame
      return(datatable)
    }
  ) %>% do.call(rbind, .)
  return(datatable)
}
