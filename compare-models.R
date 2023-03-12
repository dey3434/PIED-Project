library(loo)
library(tidyverse)
library(kableExtra)
# Load each of the models and evaluate predictive performance

MSE <- rep(NA, 10)
for (i in 1:10) {
  print(i)
  # while these models are running, ignore their values
  # if (!((i == 1) | (i == 4))) {
    # print("i=", i)
    load(here::here("results", paste0("model-", i-1, "-pred-obs.RData")))
    MSE[i] <- mean(p.o.df$error^2)
  # }
}

dat_MSE <- data.frame(model = paste0("model ", 0:9), MSE = MSE)
p_mse <- ggplot(dat_MSE, aes(x = model, y = MSE)) +
  geom_point()

ggsave(here::here("images", "compare-mse.png"), p_mse,
       width = 16, height = 9)

loos <- vector(mode = "list", length = 10)
for (i in 1:10) {
  print(i)
  # while these models are running, ignore their values
  # if (!((i == 1) | (i == 4))) {
    load(here::here("results", paste0("model-", i-1, "-loo.RData")))
    
    loos[[i]] <- leaveoneout
    rm(ll, r_eff, leaveoneout)
  # }
}

names(loos) <- paste0("model ", 0:9)
# drop the empty loos for now
# loos <- loos[lapply(loos, length) > 0]
comp <- loo::loo_compare((loos))
# comp <- loo::loo_compare(loos[[1]], loos[[2]], loos[[3]], loos[[4]], loos[[5]],
#                          loos[[6]], loos[[7]], loos[[8]])
# comp <- loo::loo_compare(loos[[1]], loos[[2]], loos[[3]], loos[[4]], loos[[5]], 
#                          loos[[6]], loos[[7]], loos[[8]], loos[[9]], loos[[10]])
print(comp, digits = 2)

print(comp, digits = 4, simplify = FALSE)

models <- rep(NA, 10)
means <- rep(NA, 10)
ses <- rep(NA, 10)
lowers <- rep(NA, 10)
uppers <- rep(NA, 10)
lowers50 <- rep(NA, 10)
uppers50 <- rep(NA, 10)
for (i in 1:10) {
  models[i]  <- paste0("model_", i-1)
  means[i] <- loos[[i]]$estimates[1, 1]
  ses[i] <- loos[[i]]$estimates[1, 2]
  lowers50[i] <- loos[[i]]$estimates[1, 1] + qnorm(0.25) * loos[[i]]$estimates[1, 2]
  uppers50[i] <- loos[[i]]$estimates[1, 1] + qnorm(0.75) * loos[[i]]$estimates[1, 2]
  lowers[i] <- loos[[i]]$estimates[1, 1] + qnorm(0.025) * loos[[i]]$estimates[1, 2]
  uppers[i] <- loos[[i]]$estimates[1, 1] + qnorm(0.975) * loos[[i]]$estimates[1, 2]
}

dat <- data.frame(model = models, mean = means, lower = lowers, upper = uppers, 
                  lower50 = lowers50, upper50 = uppers50)

p_looic <- ggplot(dat, aes(x = model, y = mean)) +
  geom_point() + 
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 1.5) +
  geom_linerange(aes(ymin = lower, ymax = upper)) +
  ylab("elpd") +
  xlab("")

ggsave(here::here("images", "compare-loo.png"), p_looic,
       width = 16, height = 9)

dat_metrics <- data.frame(
  model = models,
  MSE = MSE, 
  loo = means, 
  se = ses)

write_csv(dat_metrics, file = here::here("results", "model-results.csv"))
