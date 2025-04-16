AIC_stepwise_regression <- function(outcome, dataset, variables, stepwise_direction){
  model <- as.formula(paste0(outcome, " ~ 1 +", paste0(variables, collapse = " + ")))
  if (stepwise_direction == "forward") {
    initialize_model <- glm(formula =  as.formula(paste0(outcome, " ~ 1")), data = dataset, family = "binomial")
    stepwise_glm_model <- suppressMessages(stats::step(initialize_model, direction = stepwise_direction, scope = model, trace=0))
  }   else {
    glm_model <- glm(model, data = dataset, family = "binomial")
    stepwise_glm_model <- suppressMessages(stats::step(glm_model, direction = stepwise_direction, trace=0))
  }
  ci <- suppressMessages(confint(stepwise_glm_model))
  final <- cbind(exp(cbind(OR = coef(stepwise_glm_model), ci)) %>% round(., 2),
                 abs(summary(stepwise_glm_model)$coefficients[, "Pr(>|z|)"]) %>%
                   round(., 4)) %>% subset(rownames(.) != "(Intercept)") %>%
    `colnames<-`(c("OR", "2.5%", "97.5%", "p_value")) %>%
    as.data.frame %>% mutate(`OR (95% CI)` = paste0(OR, " (", `2.5%`, "-", `97.5%`, ")")) %>% select(`OR (95% CI)`, p_value)
  return(final)
}
