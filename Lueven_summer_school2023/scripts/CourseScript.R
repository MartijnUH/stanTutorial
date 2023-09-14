rm(list = ls())
gc()

set.seed(1954)

library(rjson)
library(bayesplot)
library(posterior)
library(ggplot2)
library(cmdstanr)
library(parallel)
library(loo)
library(outbreaks)
library(gridExtra)
source("scripts/tools_is.r")

mc.cores = detectCores()

## Linear regression model
# transpile (translate Stan to C++ and then compile)
mod <- cmdstan_model("scripts/model/linear.stan")

# run sampler
n_chains <- 4
fit <- mod$sample(data = data,
                  chains = n_chains,
                  init = init,
                  save_warmup = TRUE,
                  parallel_chains = detectCores())

# Examine Stan's default summaries
fit$summary()

# Construct diagnostic plots
pars <- c("beta", "sigma")
bayesplot::mcmc_trace(fit$draws(inc_warmup = TRUE),
                      n_warmup = 1000, pars = pars)

bayesplot::mcmc_dens_overlay(fit$draws(), pars = pars)

# Extract posterior predictive checks
yrep <- as.matrix(
  as_draws_df(fit$draws(variables = c("y_pred"))))
head(yrep)

# We don't need the chain, iteration and draw ID, so let's remove them.
yrep <- yrep[, -(11:13)]

# Plot the posterior predictions and compare it to the real data.
bayesplot::ppc_ribbon(y = data$y, yrep = yrep, x = data$x,
                      y_draw = "point") +
  theme_bw() +
  ylab("y")

## SIR Model
theme_set(theme_bw())
ggplot(data = influenza_england_1978_school) +
  geom_point(mapping = aes(x = date, y = in_bed)) +
  labs(y = "Number of students in bed")

# create a data list to be passed to Stan
cases <- influenza_england_1978_school$in_bed
N <- 763;
n_days <- length(cases)
t <- seq(0, n_days, by = 1)
t0 = 0
t <- t[-1]

#initial conditions
i0 <- 1
s0 <- N - i0
r0 <- 0
y0 = c(S = s0, I = i0, R = r0)

data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t,
                 N = N, cases = cases)

# define starting distribution
init <- function() {
  list(beta = abs(rnorm(1, mean = 2, sd = 1)),
       gamma = abs(rnorm(1, mean = 0.4, sd = 0.5)),
       phi_inv = rexp(1, rate = 5))
}

# transpile (translate Stan to C++ and then compile)
mod <- cmdstan_model("scripts/model/sir.stan")

data_sir$rel_tol <- 1e-2
data_sir$abs_tol <- 1e-2
data_sir$max_num_steps <- 1e4
n_chains <- 4
fit <- mod$sample(data = data_sir,
                  chains = n_chains,
                  parallel_chains = detectCores(),
                  init = init,
                  save_warmup = TRUE)
fit$time()

pars <- c("gamma", "beta", "phi", "R0")
fit$summary(variables = pars)


bayesplot::mcmc_trace(fit$draws(inc_warmup = TRUE),
                      n_warmup = 1000, pars = pars)
bayesplot::mcmc_dens_overlay(fit$draws(), pars = pars)

# Extract posterior predictive checks
pred_cases <- as.matrix(
  as_draws_df(fit$draws(variables = c("pred_cases"))))[, -(15:17)]

bayesplot::ppc_ribbon(y = data_sir$cases, yrep = pred_cases,
                      x = data_sir$ts, y_draw = "point") +
  theme_bw() +
  ylab("cases") + xlab("days")

# Run same model with a Poisson likelihood
mod <- cmdstan_model("scripts/model/sir_poisson.stan")

fit_poisson <- mod$sample(data = data_sir,
                          chains = n_chains,
                          parallel_chains = detectCores(),
                          init = init,
                          save_warmup = TRUE)

fit_poisson$summary(variables = pars)

pred_cases_poisson <- as.matrix(
  as_draws_df(fit_poisson$draws(variables = c("pred_cases"))))[, -(15:17)]

bayesplot::ppc_ribbon(y = data_sir$cases, yrep = pred_cases_poisson,
                      x = data_sir$ts, y_draw = "point") +
  theme_bw() +
  ylab("cases") + xlab("days")

# compute PSIS-loo estimate
log_lik_draws <- fit$draws("log_lik")
loo_estimate <- loo(log_lik_draws, r_eff = relative_eff(log_lik_draws))


log_lik_draws_poisson <- fit_poisson$draws("log_lik")
loo_estimate_poisson <-
  loo(log_lik_draws_poisson, r_eff = relative_eff(log_lik_draws_poisson))

print(loo_estimate_poisson)
print(loo_estimate)

# run SIR model with custom tolerance for the ODE solver
mod <- cmdstan_model("scripts/model/sir_tol.stan")

tol <- 1e-4
data_sir$tol <- tol
fit_tol <- mod$sample(data = data_sir,
                      chains = n_chains,
                      parallel_chains = detectCores(),
                      init = init,
                      save_warmup = TRUE)

log_ratios <- fit_tol$draws("log_ratios")

psis_fit <- psis(log_ratios, r_eff =  relative_eff(log_ratios))
psis_fit$diagnostics

# Correct Monte Carlo samplers, using importance weights.
# Only works if the log ratios don't go to 0!
is_summary(fit_tol, pars, psis_fit, log_ratios)

## 8 schools hieracichal model
n_schools <- 8
y <- c(28, 8, -3, 7, -1, 1, 18, 12)
sigma <- c(15, 10, 16, 11, 9, 11, 10, 18)

stan_data <- list(n_schools = n_schools,
                  y = y, sigma = sigma)

mod <- cmdstan_model("scripts/model/8schools_nc.stan")
# other models: 8schools_nc, 8schools_marginal

fit1 <- mod$sample(data = stan_data,
                  chains = 4, 
                  parallel_chains = detectCores(),
                  iter_warmup = 1000,
                  iter_sampling = 1000,
                  adapt_delta = 0.9,
                  seed = 1234)

pars <- c("mu", "tau")

mod <- cmdstan_model("scripts/model/8schools_marginal.stan")
# other models: 8schools_nc, 8schools_marginal

fit2 <- mod$sample(data = stan_data,
                   chains = 4, 
                   parallel_chains = detectCores(),
                   iter_warmup = 1000,
                   iter_sampling = 1000,
                   adapt_delta = 0.9,
                   seed = 1234)

# compare monte carlo estimates
pars <- c("mu", "tau")
rbind(fit1$summary(variables = pars), fit2$summary(variables = pars))

## Disease map of Finland
data <- fromJSON(file = "scripts/data/disease_100.json")
mod <- cmdstan_model("scripts/model/disease_map.stan")

fit <- mod$sample(data = data,
                  chains = 4, parallel_chains = detectCores(),
                  iter_warmup = 500,
                  iter_sampling = 500,
                  adapt_delta = 0.8,
                  seed = 1234)

pars <- c("alpha", "rho", "theta[1]")
fit$summary(variables = pars)

## Pharmacokinetic model
# Download Torsten
system("git clone https://github.com/metrumresearchgroup/Torsten.git")

# Make cmdstanr use Torsten
cmdstan_path <- "Torsten/cmdstan/"
set_cmdstan_path(cmdstan_path)
bayesplot::color_scheme_set("mix-blue-green")

# Silence certain C++ warning messages during compilation
cmdstan_make_local(
  dir = cmdstan_path(),
  cpp_options = list(
    "CXXFLAGS += -Wno-deprecated-declarations -Wno-ignored-attributes"), 
  append = TRUE)

data <- fromJSON(file = "scripts/data/twoCpt.data.json")

# Plot drug concentration profile
p <- ggplot(data = data.frame(time = data$time[-1],
                              cObs = data$cObs),
            aes(x = time, y = cObs)) +
  geom_point() +
  theme_bw()
p

# Draw initial conditions from the prior
init <- function() {
  list(CL = exp(rnorm(1, log(10), 0.25)),
       Q = exp(rnorm(1, log(15), 0.5)),
       VC = exp(rnorm(1, log(35), 0.25)),
       VP = exp(rnorm(1, log(105), 0.5)),
       ka = exp(rnorm(1, log(2), 1)),
       sigma = abs(rnorm(1, 0, 1)))
}

model_name <- "twoCpt_ode"
mod <- cmdstan_model(paste0("scripts/model/", model_name, ".stan"))

n_chains <- 4
data$tol <- 0.25
fit <- mod$sample(data = data, chains = n_chains, init = init,
                  parallel_chains = detectCores(),
                  iter_warmup = 500, iter_sampling = 500)

fit$time()
# Save fit (useful for expensive models!)
system("mkdir deliv")
fit$save_object(paste0("deliv/", model_name, ".fit.RDS"))

pars = c("CL", "Q", "VC", "VP", "ka", "sigma")
fit$summary(variables = pars)

bayesplot::mcmc_trace(fit$draws(), pars = pars)
bayesplot::mcmc_dens_overlay(fit$draws(), pars = pars)

# posterior predictive checks
yrep <- as.matrix(
  as_draws_df(
    fit$draws(variables = c("concentrationObsPred"))
  ))[, -(52:54)]

yobs <- data$cObs
time <- data$time[-1]

p <- bayesplot::ppc_ribbon(y = yobs, yrep = yrep,
                           x = time, y_draw = "point") +
  xlab("time(h)") + ylab("drug plasma concentration (mg/L)") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 18))
p

# PSIS check for ODE tolerance
fit$loo()
log_ratios <- fit$draws("log_ratios")

psis_fit <- psis(log_ratios, r_eff =  relative_eff(log_ratios))
psis_fit$diagnostics

## Use built-in analytical solution
# two compartment model
mod <- cmdstan_model("scripts/model/twoCpt.stan")
fit <- mod$sample(data = data, chains = n_chains, init = init,
                  parallel_chains = detectCores(),
                  iter_warmup = 500, iter_sampling = 500)

fit$time()
fit$summary()

yrep <- as.matrix(
  as_draws_df(
    fit$draws(variables = c("concentrationObsPred"))
  ))[, -(52:54)]

yobs <- data$cObs
time <- data$time[-1]

p <- bayesplot::ppc_ribbon(y = yobs, yrep = yrep,
                           x = time, y_draw = "point") +
  xlab("time(h)") + ylab("drug plasma concentration (mg/L)") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 18))
p

# one compartment model
mod <- cmdstan_model("scripts/model/oneCpt.stan")
fit_one <- mod$sample(data = data, chains = n_chains, init = init,
                      parallel_chains = n_chains,
                      iter_warmup = 500, iter_sampling = 500)

fit_one$time()
fit_one$summary()

yrep <- as.matrix(
  as_draws_df(
    fit_one$draws(variables = c("concentrationObsPred"))
  ))[, -(52:54)]

yobs <- data$cObs
time <- data$time[-1]

p <- bayesplot::ppc_ribbon(y = yobs, yrep = yrep,
                           x = time, y_draw = "point") +
  xlab("time(h)") + ylab("drug plasma concentration (mg/L)") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 18))
p

# estimate ELPD_loo for both models
log_lik_draws <- fit$draws("log_lik")
loo_estimate <- loo(log_lik_draws, r_eff = relative_eff(log_lik_draws))


log_lik_draws_one <- fit_one$draws("log_lik")
loo_estimate_one <-
  loo(log_lik_draws_one, r_eff = relative_eff(log_lik_draws_one))

print(loo_estimate_one)
print(loo_estimate)

## Population Pharmacokinetic Model
# plot data
patientID <- rep(NA, data$nSubjects)
for (i in 1:data$nSubjects) {
  patientID[data$start[i]:data$end[i]] <- i
}

p <- ggplot(data = data.frame(cObs = data$cObs, time = data$time[data$iObs],
                              ID = patientID[data$iObs]), aes(x = time, y = cObs)) +
  theme_bw() + geom_point() + facet_wrap(~ID)
p

# Draw initial points from the prior
init <- function () {
  n_subjects <- data$nSubjects
  pop_var <- c(0.2, 0.2, 0.2, 0.2, 0.2)
  
  CL_pop <- exp(rnorm(1, log(10), pop_var[1]))
  Q_pop <- exp(rnorm(1, log(15), pop_var[2]))
  VC_pop <- exp(rnorm(1, log(35), pop_var[3]))
  VP_pop <- exp(rnorm(1, log(105), pop_var[4]))
  ka_pop <- exp(rnorm(1, log(2.5), pop_var[5]))
  omega <- abs(rnorm(5, 0, pop_var))
  
  theta_pop <- c(CL_pop, Q_pop, VC_pop, VP_pop, ka_pop)
  theta <- matrix(NA, n_subjects, length(theta_pop))
  for (j in 1:n_subjects) {
    theta[j, ] <- exp(rnorm(length(theta_pop), log(theta_pop), omega))
  }
  
  list(CL_pop = CL_pop, Q_pop = Q_pop, VC_pop = VC_pop, VP_pop = VP_pop,
       ka_pop = ka_pop, omega = omega, theta = theta,
       sigma = abs(rnorm(1, 0, 1)))
}

model_name <- "twoCptPop"
mod <- cmdstan_model(paste0("model/", model_name, ".stan"))

n_chains <- 4
fit <- mod$sample(data = data, chains = n_chains, init = init,
                  parallel_chains = detectCores(),
                  iter_warmup = 500, iter_sampling = 500,
                  seed = 1234, adapt_delta = 0.8)

fit$save_object(paste0("deliv/", model_name, ".fit.RDS"))
fit$time()

pars = c("lp__", "CL_pop", "Q_pop", "VC_pop", "VP_pop", "ka_pop", "sigma")
fit$summary(variables = pars)
bayesplot::mcmc_trace(fit$draws(), pars = pars)
bayesplot::mcmc_dens_overlay(fit$draws(), pars = pars)

# posterior predictive checks
yrep <- as.matrix(
  as_draws_df(
    fit$draws(variables = c("concentrationObsPred"))))
yrep <- yrep[, -((ncol(yrep) - 2):ncol(yrep))]

yobs <- data$cObs
time <- data$time[data$iObs]
patientID <- with(data, rep(1:nSubjects, each = nObs / nSubjects))

# within patient predictions
bayesplot::ppc_ribbon_grouped(y = yobs, yrep = yrep, x = time, patientID,
                              y_draw = "point")

# predictions for new patient
yrepNew <- as.matrix(
  as_draws_df(
    fit$draws(variables = c("cObsNewPred"))))
yrepNew <- yrepNew[, -((ncol(yrepNew) - 2):ncol(yrepNew))]

bayesplot::ppc_ribbon_grouped(y = yobs, yrep = yrepNew, x = time, patientID,
                              y_draw = "point")

# Implementation with multiple threads per chain
model_name <- "twoCptPop_rs"
mod <- cmdstan_model(paste0("model/", model_name, ".stan"),
                     cpp_options = list(stan_threads = TRUE))

n_chains <- 1
fit_rs <- mod$sample(data = data, chains = n_chains, init = init,
                     parallel_chains = n_chains,
                     threads_per_chain = 1,
                     iter_warmup = 500, iter_sampling = 500,
                     seed = 123, adapt_delta = 0.8)

fit_rs$time()
# fit_rs$summary(variables = pars)