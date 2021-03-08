# code for stan portion of the model
# must run Read_PIED_data_format.R before

sink("pied_grow.stan")
cat("
    data {
    
    int<lower=0> K;         // N. predictors 
    int<lower=0> nG;        // N. observations
    int<lower=0> nGtest;    // N. observations (ppc)
    matrix[nG,K] xG;        // Predictor matrix
    matrix[nGtest, K] xGtest;   // Predictor matrix (ppc)
    vector[nG] yG;          // log size at time t+1 
    
    int<lower=0> nplot;         // number of plots
    int<lower=1> plot[nG];      // index for plot
    
    int<lower=0> ntree;          // number of trees
    int<lower=1> tree[nG];          // index for trees
    int<lower=1> treetest[nGtest];  //index for trees (ppc)
    int<lower=1> plotfortree[ntree];   // plot index for each tree
    }
    
    parameters {
    
    real u_beta0;                          // intercept means
    vector[K] u_beta;                      // other coeff mean
    
    real beta0_p_tilde[nplot];                   // plot-level intercepts
    real<lower=0> s_beta0_p;               // plot-level intercept variance
    real beta0_t_tilde[ntree];                   // tree-level intercepts
    real<lower=0> s_beta0_t;               // tree-level intercept variance
    
    real<lower=0> sigma_y;                 // Residual for growth model
    
    }
    
    transformed parameters {
    real beta0_p[nplot];                   // plot-level intercepts
    real beta0_t[ntree];                   // tree-level intercepts
    
    for(p in 1:nplot){
    beta0_p[p] = u_beta0 + s_beta0_p * beta0_p_tilde[p];
    }
    
    for(t in 1:ntree){
    beta0_t[t] = beta0_p[plotfortree[t]] + s_beta0_t * beta0_t_tilde[t];
    }
    
    }
    
    model {
    vector[nG] mG;
    
    u_beta0 ~ normal(0, 100);
    beta0_p_tilde ~ normal(0,1);
    beta0_t_tilde ~ normal(0,1);
    
    u_beta ~ normal(0, 100); 
    s_beta0_p ~ cauchy(0,2.5);
    s_beta0_t ~ cauchy(0,2.5);
    
    sigma_y ~ gamma(2,0.01);
    
    //tried nesting random effects of trees within plots--had issue
    //for(p in 1:nplot){
    //beta0_p[p] ~ normal(u_beta0, s_beta0_p);
    //for(t in 1:ntree){
    //beta0_t[t] ~ normal(beta0_p[plotfortree[t]], s_beta0_t);
    //}
    //}
    
    
    // GROWTH MODEL
    
    for(n in 1:nG){
    mG[n] = beta0_t[tree[n]]+xG[n]*u_beta;
    }
    
    yG ~ normal(mG,sigma_y);
    //yG ~ gamma(mG,sigma_y);

    }
    
    generated quantities{
    vector[nGtest] yrep;
    for(n in 1:nGtest){
    yrep[n] = normal_rng(beta0_t[treetest[n]]+xGtest[n]*u_beta,sigma_y);
    }

    //for(n in 1:nGtest){
    //yrep[n] = gamma_rng(beta0_t[treetest[n]]+xGtest[n]*u_beta,sigma_y);
    //}

    }
    
    ",fill=T)

sink()


fit_grow <- stan(file = 'pied_grow.stan', data = pied_dat, 
                 iter = 5000, warmup = 1000, chains = 3)
#saveRDS(fit_grow, file = "log_normal_fg_7_24_20.RDS")
saveRDS(fit_grow, file = "ppt_tmp_springfall_sizefix.RDS")
summary<-summary(fit_grow)
summary

plotdata<-select(as.data.frame(fit_grow),"yrep[1]":"yrep[8780]")
plotdatainterval<-select(as.data.frame(fit_grow), "u_beta[1]":"u_beta[28]")
colnames(plotdatainterval) <- c("u_beta_ppt_norm", "u_beta_tmp_norm", "u_beta_Precip_JulAug", "u_beta_Precip_NovDecJanFebMar", 
                                "u_beta_Tmean_AprMayJun", "u_beta_Tmean_SepOct", "u_beta_DIA_prev",
                                "u_beta_ppt_norm_tmp_norm", "u_beta_ppt_norm_Precip_JulAug", 
                                "u_beta_ppt_norm_DIA_prev", "u_beta_ppt_norm_Precip_NovDecJanFebMar",
                                "u_beta_ppt_norm_Tmean_AprMayJun", "u_beta_ppt_norm_Tmean_SepOct",
                                "u_beta_tmp_norm_DIA_prev", "u_beta_tmp_norm_Precip_JulAug",
                                "u_beta_tmp_norm_Precip_NovDecJanFebMar", "u_beta_tmp_norm_Tmean_AprMayJun",
                                "u_beta_tmp_norm_Tmean_SepOct", "u_beta_DIA_prev_Precip_JulAug",
                                "u_beta_DIA_prev_Precip_NovDecJanFebMar", "u_beta_DIA_prev_Tmean_AprMayJun",
                                "u_beta_DIA_prev_Tmean_SepOct", "u_beta_Precip_JulAug_Precip_NovDecJanFebMar",
                                "u_beta_Precip_JulAug_Tmean_AprMayJun", "u_beta_Precip_JulAug_Tmean_SepOct",
                                "u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun", "u_beta_Precip_NovDecJanFebMar_Tmean_SepOct",
                                "u_beta_Tmean_AprMayJun_Tmean_SepOct")
ppc_dens_overlay(yGtest, as.matrix(plotdata))

ext_fit <- rstan::extract(fit_grow)
yrep <- ext_fit$yrep
#yrep <- exp(yrep)
mean.pred <- apply(ext_fit$yrep, 2, median)
p.o.df <- data.frame(predicted = exp(mean.pred), observed = exp(grow_test$loggrowth), error = (exp(mean.pred) - exp(grow_test$loggrowth)))
meansqrd <- (mean(p.o.df$error))^2
ggplot(p.o.df, aes(predicted, observed)) + geom_point(alpha = 0.1) + geom_abline(aes(intercept = 0, slope = 1), color = "red", linetype = "dotted") +
  ylim(0, 10) + xlim(0,10)


#Validation
sigma <- as.data.frame(fit_grow)[,"sigma_y"]
mu <- as.matrix(plotdatainterval) %*% t(xG)
ll <- matrix(0, length(sigma), length(yG))
for(i in 1:length(sigma)){
  ll[i,] <- dnorm(yG, mu[i,], sd = sigma[i], log = TRUE)
}
newll <- as.matrix(ll)
r_eff <- relative_eff(exp(ll), chain_id = rep(1:3, each = 4000), cores = 8)
leaveoneout <- loo(as.matrix(ll), r_eff = r_eff, save_psis = TRUE, cores = 8)

yrep <- matrix(0, length(sigma), length(yG))
for(i in 1:length(sigma)){
  yrep[i,] <- rnorm(yG, mu[i,], sd = sigma[i])
}
psis <- leaveoneout$psis_object
keep_obs <- sample(1:length(yG), 100)
lw <- weights(psis)
ppc_loo_intervals(yG, yrep = yrep, psis_object = psis, subset = keep_obs, order = "median") 
ppc_loo_pit_overlay(yG, yrep = yrep, lw = lw)
ppc_loo_pit_qq(yG, yrep = yrep, lw = lw)


yrep_mean <- apply(yrep, MARGIN = 2, FUN = mean)
yrep_ci.low <- apply(yrep, MARGIN = 2, FUN = function(x){quantile(x, 0.025)})
yrep_ci.high <- apply(yrep, MARGIN = 2, FUN = function(x){quantile(x, 0.975)})

p.o.df <- data.frame(ci.low = yrep_ci.low, mean = yrep_mean, ci.high = yrep_ci.high, observed = yG,
                     ppt_norm = grow_train$ppt_norm, tmp_norm = grow_train$tmp_norm, Precip_JulAug = grow_train$Precip_JulAug, 
                     Precip_NovDecJanFebMar = grow_train$Precip_NovDecJanFebMar, Tmean_AprMayJun = grow_train$Tmean_AprMayJun, 
                     Tmean_SepOct = grow_train$Tmean_SepOct, DIA_prev = grow_train$DIA_prev, ELEV = grow_train$ELEV, 
                     SLOPE = grow_train$SLOPE, ASPECT = grow_train$ASPECT, tmp_yr = grow_train$tmp_yr, ppt_yr = grow_train$ppt_yr, 
                     year = grow_train$year, PlotCD = grow_train$PlotCD, treeCD = grow_train$treeCD,
                     name = grow_train$name)
p.o.df$overpredicted <- ifelse(p.o.df$observed <= p.o.df$ci.low, "over predicted", "within confidence interval")
overpredictedpoints <- p.o.df %>% filter(observed <= ci.low)

# ggplot(data = p.o.df, aes(x = ELEV, y = observed, color = overpredicted)) + geom_point(size = 0.5)
# ggplot(data = p.o.df, aes(x = tmp_norm, y = observed, color = overpredicted)) + geom_point(size = 0.5)
# ggplot(data = p.o.df, aes(x = ppt_norm, y = observed, color = overpredicted)) + geom_point(size = 0.5)
# ggplot(data = p.o.df, aes(x = Precip_JulAug, y = observed, color = overpredicted)) + geom_point(size = 0.5)
# ggplot(data = p.o.df, aes(x = Precip_NovDecJanFebMar, y = observed, color = overpredicted)) + geom_point(size = 0.5)
# ggplot(data = p.o.df, aes(x = Tmean_AprMayJun, y = observed, color = overpredicted)) + geom_point(size = 0.5)
# ggplot(data = p.o.df, aes(x = Tmean_SepOct, y = observed, color = overpredicted)) + geom_point(size = 0.5)
# ggplot(data = p.o.df, aes(x = DIA_prev, y = observed, color = overpredicted)) + geom_point(size = 0.5)
# ggplot(data = p.o.df, aes(x = SLOPE, y = observed, color = overpredicted)) + geom_point(size = 0.5)
# ggplot(data = p.o.df, aes(x = ASPECT, y = observed, color = overpredicted)) + geom_point(size = 0.5)
# ggplot(data = p.o.df, aes(x = tmp_yr, y = observed, color = overpredicted)) + geom_point(size = 0.5)
# ggplot(data = p.o.df, aes(x = ppt_yr, y = observed, color = overpredicted)) + geom_point(size = 0.5)
# ggplot(data = p.o.df, aes(x = year, y = observed, color = overpredicted)) + geom_point(size = 0.5)
# ggplot(data = p.o.df, aes(x = name, y = observed, color = overpredicted)) + geom_point(size = 0.5)

#----------------------------------------------------------------------------------------------
#  Run other validation code!
#----------------------------------------------------------------------------------------------

