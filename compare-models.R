# Load each of the models and evaluate predictive performance

MSE <- rep(0, 9)
for (i in 1:9) {
  load(here::here("results", paste0("model-", i, "-pred-obs.RData")))
  MSE[i] <- mean(p.o.df$error^2)
}


loos <- vector(mode = "list", length = 9)
for (i in 1:9) {
  load(here::here("results", paste0("model-", i, "-loo.RData")))
}