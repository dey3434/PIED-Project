# Load each of the models and evaluate predictive performance

MSE <- rep(NA, 10)
for (i in 1:10) {
  print(i)
  # while these models are running, ignore their values
  if (!((i == 1) | (i == 4))) {
    print("i=", i)
    load(here::here("results", paste0("model-", i-1, "-pred-obs.RData")))
    MSE[i] <- mean(p.o.df$error^2)
  }
}

dat <- data.frame(model = paste0("model_", 0:9), MSE = MSE)
ggplot(dat, aes(x = model, y = MSE)) +
  geom_point()

loos <- vector(mode = "list", length = 10)
for (i in 1:10) {
  print(i)
  # while these models are running, ignore their values
  if (!((i == 1) | (i == 4))) {
    load(here::here("results", paste0("model-", i-1, "-loo.RData")))
    
    loos[[i]] <- leaveoneout
    rm(ll, r_eff, leaveoneout)
  }
}


# drop the empty loos for now
loos <- loos[lapply(loos, length) > 0]
comp <- loo::loo_compare((loos))
# comp <- loo::loo_compare(loos[[1]], loos[[2]], loos[[3]], loos[[4]], loos[[5]],
#                          loos[[6]], loos[[7]], loos[[8]])
# comp <- loo::loo_compare(loos[[1]], loos[[2]], loos[[3]], loos[[4]], loos[[5]], 
#                          loos[[6]], loos[[7]], loos[[8]], loos[[9]], loos[[10]])
print(comp, digits = 2)

print(comp, digits = 4, simplify = FALSE)
