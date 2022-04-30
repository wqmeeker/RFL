
###
###  stand-alone RFL github code
###

library(rstan)
options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###
### read the data and make the Stan data list
###   including highly informative priors around the ML estimates
###
LaminatePanel <-  read.csv("LaminatePanel.csv")
head(LaminatePanel)

LaminatePanelRFL.dat <- list(
N = nrow(LaminatePanel),
M = 1,
tran_response = as.matrix(log(LaminatePanel[, "Thousands.of.Cycles"])),
cencode = ifelse(as.character(LaminatePanel[, "Censoring.Indicator"])=="Failed", 1, 2),
weights = rep(1.0, length=nrow(LaminatePanel)),
stress = LaminatePanel[, "Stress..MPa."],
log_s0_minus_gamma0 = 0,
numdist = 4,
prior_codes = c(1,1,1,1,1),
beta0star_prior_parameters = c(29.,  2.329, 0),
beta1_prior_parameters = c(-4.75, 0.0970, 0),
log_sigma_error_prior_parameters = c(-0.948, 0.099, 0),
mu_log_gamma_prior_parameters = c(1.70, 0.02828445, 0),
log_sigma_log_gamma_prior_parameters = c(-4.71, 0.312, 0))
 
###
### compile the Stan RFL model
###                         )
RFL <- stan_model(file='RFL.stan')

###
### get the initial values list
###
source("sim.init.R")
source("rangePrior.R")
the.inits <- sim.init(param.names=c("beta0star", "beta1", "log_sigma_error", "mu_log_gamma", "log_sigma_log_gamma"), param.lower=c(27., -4.6, log(0.46), 5.3, log(0.009)), param.upper=c(29., -4.5, log(0.47), 5.4, log(0.01)), distribution=c("normal", "normal", "normal", "normal", "normal"), number=4)

###
###  run the estimation
###
LaminatePanelRFL.fit <- sampling(RFL, data = LaminatePanelRFL.dat, iter = 10000, chains = 4, init=the.inits, control = list(max_treedepth = 14, adapt_delta=0.99), thin=5, warmup=2000, open_progress=TRUE)

###
### inspect the ressults
###
the.parameters <-  c("beta0star", "beta1", "log_sigma_error", "mu_log_gamma", "log_sigma_log_gamma")
interesting.pars <- c("beta0", "beta1", "sigma_error", "mu_log_gamma", "sigma_log_gamma")

print(LaminatePanelRFL.fit, pars=the.parameters, digits=3)
pairs(LaminatePanelRFL.fit, pars=the.parameters)
traceplot(LaminatePanelRFL.fit, pars=the.parameters)

print(LaminatePanelRFL.fit, pars=interesting.pars, digits=3)
pairs(LaminatePanelRFL.fit, pars=interesting.pars)
traceplot(LaminatePanelRFL.fit, pars=interesting.pars)

