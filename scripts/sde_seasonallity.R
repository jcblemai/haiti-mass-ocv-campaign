# Title: Seasonality of diarrheal cases in BF
# Description: Exploring seasonal patterns of diarrhea using SDE filtering
# Date: Tue Oct 23 08:05:25 2018
# Author: javier.perezsaez@epfl.ch



# Preamble ---------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(pomp)
library(foreach)
library(doSNOW)
library(magrittr)
library(lubridate)
rm(list = ls())
Sys.setlocale("LC_ALL","C")

# function to convert dates to fractions of years for model
dateToYears <- function(date, origin = as.Date("2014-01-01"), yr_offset = 2014) {
  julian(date, origin = origin)/365.25 + yr_offset
}

yearsToDate <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
  as.Date((year_frac - yr_offset) * 365.25, origin = origin)
}

yearsToDateTime <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
  as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
}

# function to draw random vriables from a Drichillet distribution using independent gammas
rdirichlet <- function (n, alpha) {
  if (!is.null(salpha <- dim(alpha))) {
    x <- matrix(rgamma(prod(salpha) * n, rep(as.vector(as.matrix(alpha)), each = n)), ncol = salpha[2], byrow = F)
    sm <- x %*% rep(1, salpha[2])
  } else {
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
  }
  return(x/as.vector(sm))
}

# Load data ---------------------------------------------------------------

# Test for Artibonite
department <- "Artibonite"

# Load the weekly cholera cases
cases <- read_csv("haiti-data/fromAzman/cases_corrected.csv")  %>% 
  gather(dep, cases, -date) %>% 
  filter(dep == department) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date))

# get the time of the first datapoint (use %>% filter(time > 2015) to constraint)
t_first_datapnt <- cases %>% slice(1) %>% .[["time"]]

# Load the rainfall
rain <- read_csv("haiti-data/fromAzman/rainfall.csv")  %>% 
  gather(dep, rain, -date) %>% 
  group_by(dep) %>% 
  mutate(max_rain = max(rain), rain_std = rain/max_rain) %>%
  ungroup() %>% 
  filter(dep == department) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date))


# plot data
cases %>%
  select(date, cases) %>% 
  gather(variable, value, -date) %>% 
  bind_rows(
    rain %>%
      select(date, rain) %>% 
      gather(variable, value, -date)
  ) %>% 
  mutate(yr = year(date), month = month(date)) %>% 
  filter(yr > 2014, yr < 2016, month > 4, month < 7) %>% 
  ggplot(aes(x = date, y = value)) +
  geom_line() +
  facet_grid(variable~., scales = "free_y")

### Model specification -----------------------------------------------------
## state variable names:
# S:  Susceptibles
# I:  Infected
# R:  Recovered + asymptomatic infected

## data names:
# cases: reported suspected cholera cases (monthly)

## paramter names:
### Pop dynamics
# H:        total population
# gamma:    recovery from infection
# mu:       natural mortylity
# phi:      fertility
# alpha:    diarrhea-induced mortality
# rho:      loss of aquired immunity

### Infection dynamics

#### Basic SIRB
# beta:   directed (human-to-human) transmission coefficient

#### Extra-demographic stochasticity
# std_W:    standard deviation of the weiner process to perturb the foi       

### Measurement model
# epsilon:  under-rerpoting fraction
# k: over-dispertion paramter of negative binomial

# Set variables -----------------------------------------------------------

# define stat variable names
state_names <- c("S", "I", "R", "C", "X")

# define parameter names for pomp
## process model parameters names to estimate
param_proc_est_names <- c("rho", "sigma_X", "epsilon", "k", "std_W")

## initial value parameters to estimate
param_iv_est_names <- c("S_0", "I_0", "R_0", "X_0")

## fixed process model parameters 
param_proc_fixed_names <- c("H",  "mu", "phi", "gamma", "alpha")

## fixed initial value parameters
param_iv_fixed_names <- NULL

# all paramter names to estimate
param_est_names <- c(param_proc_est_names, param_iv_est_names)
# all fixed parameters
param_fixed_names <- c(param_proc_fixed_names, param_iv_fixed_names)

# all param names
param_names <- c(param_est_names, param_fixed_names)

# names of parameters that are rates
param_rates_in_days_names <- c("mu", "phi", "gamma", "rho")

# names of rate parameters
param_rate_names <- param_names[!str_detect(param_names, "_0|H|k|epsilon|eff|t_|std_W|sigma|alpha_|lambda")]


# Measurment model  -------------------------------------------------------

# measurement model
## density

## NegBinomial density (if k -> inf then becomes Poisson)
dmeas <- Csnippet("double mean_cases;
                  mean_cases = epsilon * C;
                  
                  if (ISNA(cases)) {
                    lik = (give_log) ? 0 : 1;
                  } else {
                    lik = dnbinom_mu(cases, k, mean_cases, give_log) ;
                  }
                  ")

## NegBinomial simulator
rmeas <- Csnippet("
                  double mean_cases = epsilon * C;
                  cases = rnbinom_mu(k, mean_cases);
                  ")

# Process model -----------------------------------------------------------------


sirb.rproc <- Csnippet("
                       double foi, foi_stoc; // force of infection and its stochastic version
                       double dw, dWX;       // extra-demographic stochasticity on foi
                       double rate[6];      // vector of all rates in model
                       double dN[6];        // vector of transitions between classes during integration timestep
                       int births;    // number of births in time step

                       // force of infection
                       foi = exp(X) * I / (S+I+R);
                       
                       if(std_W > 0.0) {
                       // white noise (extra-demographic stochasticity)
                       dw = rgammawn(std_W, dt);
                       // apply stochasticity
                       foi_stoc = foi * dw/dt;
                       } else {
                       foi_stoc = foi;
                       }
                       
                       // define transition rates for each type of event
                       // S compartment
                       rate[0] = foi_stoc;   // infections
                       rate[1] = mu;         // natural deaths
                       // I compartment
                       rate[2] = mu;         // natural deaths
                       rate[3] = alpha;      // cholera-induced deaths
                       rate[4] = gamma;      // recovery from infection
                       // R compartment
                       rate[5] = rho;        // loss of natural immunity
                       rate[6] = mu;         // natural death
                      
                       births = rpois((S+I+R) * phi * dt);
                       // simulate all transitions
                       reulermultinom(2, S, &rate[0], dt, &dN[0]);
                       reulermultinom(3, I, &rate[2], dt, &dN[2]);
                       reulermultinom(2, R, &rate[5], dt, &dN[5]);
                       
                       // update state variables
                       //S   += births - dN[0] - dN[1] + dN[5];
                       I   += dN[0] - dN[2] - dN[3] - dN[4];
                       R   += dN[4] - dN[5] - dN[6];
                       C   += dN[0];

                       S = nearbyint(H - I - R);
                       // random walk of beta
                       dWX = rnorm(0, sqrt(dt));
                       X  +=  sigma_X * dWX;
                       ")


# Initializer -------------------------------------------------------------

initalizeStates <- Csnippet("
                            double m = H /(S_0+I_0+R_0);
                            S   = nearbyint(S_0 * m);
                            I   = nearbyint(I_0 * m);
                            R   = nearbyint(R_0 * m);
                            C   = 0;
                            X   = X_0;
                            ")


# Parameter transformations -----------------------------------------------
# use log for positive parameters and logit for parmaters in [0,1]

toEstimationScale <- Csnippet("
                              Tsigma_X = log(sigma_X);
                              Tstd_W = log(std_W);
                              Tgamma = log(gamma);
                              Talpha = log(alpha);
                              Trho = log(rho);
                              Tepsilon = logit(epsilon);
                              Tk = log(k);
                              TS_0 = logit(S_0);
                              TI_0 = logit(I_0);
                              TR_0 = logit(R_0);
                              ")

fromEstimationScale <- Csnippet("
                                Tsigma_X = exp(sigma_X);
                                Tstd_W = exp(std_W);
                                Tgamma = exp(gamma);
                                Talpha = exp(alpha);
                                Trho = exp(rho);
                                Tepsilon = expit(epsilon);
                                Tk = exp(k);
                                TS_0 = expit(S_0);
                                TI_0 = expit(I_0);
                                TR_0 = expit(R_0);
                                ")


# Build pomp object -------------------------------------------------------

# input parameters to the model
input_parameters <- yaml::read_yaml("haiti-data/input_parameters_sde.yaml")

# Start and end dates of epidemic
t_start <- dateToYears(as.Date(input_parameters$t_start))
t_end <- dateToYears(as.Date(input_parameters$t_end))

# get fixed process paramteres to input
fixed_input_parameters <- as_vector(unlist(input_parameters[map_lgl(names(input_parameters), ~ . %in% param_fixed_names)]))

# set fixed process parameters
param_proc_fixed <- set_names(seq_along(param_proc_fixed_names) * 0, param_proc_fixed_names)
param_proc_fixed[names(fixed_input_parameters)] <- fixed_input_parameters

# Initialize the fixed parameters
param_fixed <-  set_names(seq_along(param_fixed_names) * 0, param_fixed_names)
param_fixed[param_proc_fixed_names] <- as.numeric(param_proc_fixed)

# The population parameter depends on the health center, needs to be initialized when fitting
pops <- yaml::read_yaml("haiti-data/input_parameters.yaml")$population  # population of all departments
param_fixed["H"] <- unlist(pops)[department]


# Initialize the parameters to estimate (just initial guesses)
param_est <- set_names(seq_along(param_est_names) * 0, param_est_names)
param_est["rho"] <- 1/(365 * 3)
param_est["gamma"] <- 1/(5)
param_est["alpha"] <- 0
param_est["std_W"] <- 0
param_est["epsilon"] <- .2
param_est["k"] <- 10
param_est["R_0"] <- 0.1
param_est["S_0"] <- 0.8
param_est["I_0"] <- 0.1
param_est["X_0"] <- 3.4
param_est["sigma_X"] <- 1

# rate of simulation in fractions of years
dt_yrs <- 1/365.25

# adjust the rate parameters depending on the integration delta time in years (some parameter inputs given in days)
params <- c(param_est, param_fixed)
params[param_rates_in_days_names] <- params[param_rates_in_days_names] * 365.25

# Setup MIF paramters -----------------------------------------------------

# values of the random walks standard deviations
rw.sd_rp <- 0.02  # for the regular (process and measurement model) paramters
rw.sd_ivp <- 0.2  # for the initial value paramters
rw.sd_param <- set_names(c(rw.sd_rp, rw.sd_ivp), c("regular", "ivp"))

# Level of detail on which to run the computations
run_level <- 2
sir_Np <-           c(1000, 3e3, 1e4)
sir_Nmif <-         c(1,    200,  400)
sir_Ninit_param <-  c(8,    8,  10)
sir_NpLL <-         c(1000, 1e4, 5e4)
sir_Nreps_global <- c(1,    10,   100)

# lower bound for positive parameter values
min_param_val <- 1e-5 
# define the bounds for the paramters to estimate
parameter_bounds <- tribble(
  ~param, ~lower, ~upper,
  # Regular paramters
  # "gamma", min_param_val, 3,
  # "alpha", min_param_val, 2,
  "rho", min_param_val, 100,
  # Process noise
  "sigma_X", min_param_val, 5,
  # "std_W", min_param_val, 1e-1,
  # Measurement model
  "epsilon", min_param_val, 1,
  "k", min_param_val, 10,
  # Initial conditions
  "S_0", min_param_val, 1,
  "I_0", min_param_val, .1,
  "R_0", min_param_val, .1,
  "X_0", min_param_val, 3
)

# convert to matrix for ease
parameter_bounds <- set_rownames(as.matrix(parameter_bounds[, -1]), parameter_bounds[["param"]])

# create random vectors of initial paramters given the bounds
init_params <- sobolDesign(lower = parameter_bounds[, "lower"],
                           upper = parameter_bounds[, "upper"], 
                           nseq =  sir_Ninit_param[run_level]) 

# set the random values of the initial conditions paramters so as they sum to 1, given the sobol design values of S, I and R
init_params[, c("S_0", "I_0", "R_0")] <- rdirichlet(1, init_params[, c("S_0", "I_0", "R_0")])

# bind with the fixed valued paramters
init_params <- cbind(init_params, 
                     matrix(rep(params[param_fixed_names],
                                each = sir_Ninit_param[run_level]),
                            nrow = sir_Ninit_param[run_level]) %>% 
                       set_colnames(param_fixed_names))

job_rw.sd <- eval(
  parse(
    text = str_c("rw.sd(",
                 # "gamma  = ",  rw.sd_param["regular"],
                 # ", alpha  = ",   rw.sd_param["regular"],
                 " rho  = ",  rw.sd_param["regular"],
                 ", sigma_X  = ",  rw.sd_param["regular"],
                 # ", std_W  = ", rw.sd_param["regular"],
                 ", epsilon  = ",  rw.sd_param["regular"],
                 ", k  = ",  rw.sd_param["regular"],
                 ", S_0  = ivp(",  rw.sd_param["ivp"], ")",
                 ", I_0  = ivp(",  rw.sd_param["ivp"], ")",
                 ", R_0  = ivp(",  rw.sd_param["ivp"], ")",
                 ", X_0  = ivp(",  rw.sd_param["ivp"], ")",
                 ")")
  )
)

# files for results
ll_filename <- "results/loglik_sde_seasonality.csv"

# parallel computations
cl <- makeCluster(parallel::detectCores())
registerDoSNOW(cl)

data <- cases %>%
  filter(time > t_start) %>% 
  select(time, cases)

sir_sde <- pomp(
  # set data
  data =  data, 
  # time column
  times = "time",
  # initialization time
  t0 = t_start - dt_yrs,
  # paramter vector
  params = params,
  # process simulator
  rprocess = euler.sim(step.fun = sirb.rproc, delta.t = dt_yrs),
  # measurement model simulator
  rmeasure =  rmeas,
  # measurement model density
  dmeasure = dmeas,
  # names of state variables
  statenames = state_names,
  # names of accumulator variables to be re-initalized at each observation timestep 
  # (C for cases)
  zeronames = c("C"),
  # names of paramters
  paramnames = param_names,
  initializer = initalizeStates,
  toEstimationScale = toEstimationScale,
  fromEstimationScale = fromEstimationScale
)

save()

# MIF it! 
t1 <- system.time({
  mf <- foreach(parstart = iter(init_params, by = "row"),
                .inorder = T, 
                .packages = "pomp") %dopar% {
                  
                  mif2(sir_sde,
                       start = parstart,
                       Np = sir_Np[run_level],
                       Nmif = sir_Nmif[run_level],
                       cooling.type = "geometric",
                       cooling.fraction.50 = 0.9,
                       transform = TRUE,
                       rw.sd = job_rw.sd,
                       verbose = F
                  )
                  
                }
})
save(mf, file = sprintf("results/mif_tests_sir_sde_%s.rda", department))

# Filter ------------------------------------------------------------------
# sir_filtered <- pfilter(sir_sde, coef(mf), Np = 10000, save.states = T)

t2 <-  system.time({
  liks <- foreach(mfit = mf,
                  .inorder = T, 
                  .packages = "pomp",
                  .combine = rbind) %dopar% {
                    
                    # compute log-likelihood estimate by repeating filtering with the given param vect
                    ll <-  logmeanexp(
                      replicate(sir_Nreps_global[run_level],
                                logLik(
                                  pfilter(sir_sde,
                                          params = pomp::coef(mfit),
                                          Np = sir_NpLL[run_level])
                                )
                      ), se = TRUE)
                    
                    # save to dataframe
                    data.frame(loglik = ll[1], loglik_se = ll[2], t(coef(mfit)))
                  } %>% 
    mutate(department = department)
})

write_csv(liks, path = ll_filename, append = file.exists(ll_filename))


# stop cluster
stopCluster(cl)
closeAllConnections()


# Filter states -----------------------------------------------------------
# 
# 
# # Run filtering
param_names <- names(coef(sir_sde))
best_loglik <- read_csv(ll_filename, col_names = T) 

best_loglik %>% 
  arrange(desc(loglik))
# sort logliks
sorted_liks <- sort(best_loglik$loglik, decreasing = T)



# Perform filtering of time series
best_param <- best_loglik[best_loglik$loglik == sorted_liks[3], param_names]
# best_param["std_W"] <- 0
sir_filtered <- pfilter(sir_sde, best_param, 1e4, save.states = T, filter.mean = T)

# recover quantiles of states
filtered <- foreach(tstep = sir_filtered@saved.states, t = sir_filtered@times, .combine = rbind) %do% {
  df <- map_df(
    rownames(tstep),
    ~data.frame(time = t, var = ., t(quantile(tstep[.,], c(0.025, 0.25, 0.5, 0.75, 0.975))))
  ) %>%
    set_colnames(c("time", "var", "q025", "q25", "q50", "q75", "q975"))
  # Take the exp of beta
  df[df$var == "X", -c(1:2)] <- exp(df[df$var == "X", -c(1:2)])
  df
} %>%
  as_tibble() %>%
  mutate(date = yearsToDate(time),
         department = department)

# write_csv(filtered, path = "results/filtered_states_2.csv", append = file.exists("results/filtered_states_2.csv"))

# add values of mean filtered cases
filtered %<>%
  bind_rows(filtered %>%
              filter(var == "C") %>%
              gather(qnt, val, -time, -var, -date, -department) %>%
              mutate(cases = coef(sir_filtered)[["epsilon"]] * val)  %>%
              select(time, date, qnt, cases) %>%
              spread(qnt, cases) %>%
              mutate(var = "cases")
  )

# plot filtered random walk of beta
p <- ggplot(filtered, aes(x = date)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.1) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3) +
  geom_line(aes(y = q50)) +
  geom_point(data = tibble(date = yearsToDate(sir_filtered@times, origin = as.Date("2013-01-01"), yr_offset = 2013),
                           cases = as.vector(sir_filtered@data),
                           var = "cases"),
             aes(x = date, y = cases), pch = 1) +
  facet_wrap(~var, scales = "free_y")
p

rain %>% 
  mutate(wk = week(date),
         yr = year(date)) %>%
  # group_by(yr, wk) %>% 
  # summarize(rain = sum(val),
  #           max_rain = max(val), 
  #           n = n()) %>% 
  # ungroup %>% 
  # filter(yr > 2014) %>% 
  right_join(filtered %>% 
               filter(var == "X") %>% 
               mutate(wk = week(date),
                      yr = year(date)),
             by = c("yr", "wk")) %>% 
  filter(time.x == time.y) %>% 
  select(date.x, q50, rain) %>% 
  gather(variable, value, -date.x) %>% 
  ggplot(aes(x = date.x, y = value)) +
  geom_line()+
  facet_grid(variable~., scales = "free_y")



filtered %>% 
  filter(var == "X" | var == "S") %>% 
  gather(quantile, value, -time, -var, -date, -department) %>% 
  spread(var, value) %>% 
  filter(quantile == "q50") %>% 
  ggplot(aes(x = date, y = S * X/coef(sir_sde)["H"]/(coef(sir_sde)["mu"]+coef(sir_sde)["gamma"]))) +
  geom_line() +
  geom_hline(aes(yintercept = 1))


filtered %>% 
  filter(var == "X") %>% 
  mutate(wk = week(date),
         yr = year(date-120)) %>% 
  group_by(yr) %>% 
  mutate(julian = julian(date, as.Date(str_c(min(yr), 1, 1, sep = "-")))
         , q50 = q50 - max(q50)
         ,julian = julian - julian[which.max(q50)]
         ) %>% 
  ggplot(aes(x = julian, y = q50, color = factor(yr))) +
  geom_line()


# ggsave(p, filename = sprintf("results/figures/filtered_states_%s.pdf", id), width = 10, height = 6)

# # plot the seasonality of beta
# filtered_all <- foreach(tstep = sir_filtered@saved.states, t = sir_filtered@times, .combine = rbind) %do% {
#   data.frame(time = t,  value = exp(tstep[5,]))
# } %>%
#   as_tibble() %>%
#   mutate(date = yearsToDate(time),
#          year = year(date),
#          department = id)
# 
# # write_csv(filtered_all, path = "results/filtered_all_beta.csv", append = T)
# 
# p2 <- ggplot(filtered_all, aes(x = month(date), y = value)) +
#   geom_smooth(aes(color = factor(year), fill = factor(year)))
# 
# ggsave(p2, filename = sprintf("results/figures/beta_seasonality_%s.pdf", id), width = 8, height = 6)
# # Simulate ----------------------------------------------------------------
# nsim <- 1
# seed <- 101
# # coef(sir_sde) <- coef(mf)
# coef(sir_sde) <- coef(a)
# pomp::simulate(sir_sde, nsim = nsim, as.data.frame = T , include.data = F, seed = seed) %>%
#   as_tibble() %>%
#   mutate(isdata = sim == "data") %>%
#   gather(variable, value, -time, -sim, -isdata) %>%
#   group_by(time, isdata, variable) %>%
#   summarise( q05 = quantile(value, 0.025, na.rm = T),
#              mean = mean(value, na.rm = T),
#              q50 = quantile(value, 0.5, na.rm = T),
#              q95 = quantile(value, 0.975, na.rm = T)) %>%
#   ungroup %>%
#   mutate(isdata = ifelse(isdata, "data", "simulation"),
#          date = yearsToDateTime(time)) %>%
#   ggplot(aes(x = date, color = isdata, fill = isdata)) +
#   geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.1) +
#   geom_line(aes(y = mean)) +
#   facet_wrap(~variable, scales = "free_y")
# # 
# 
# states <- read_csv("results/filtered_states.csv", col_names = F) %>% 
#   set_colnames(c("time", "var", "q025", "q25", "q50", "q75", "q975", "date", "department")) %>% 
#   filter(var == "X") %>% 
#   gather(qnt, value, -time, -var, -date, -department) %>% 
#   group_by(department, qnt) %>% 
#   mutate(mval = mean(value),
#          value = value/mval) %>% 
#   select(-mval) %>% 
#   spread(qnt, value)
# 
# ggplot(states, aes(x = month(date), fill = department, color = department)) +
#   geom_smooth(aes(y = q50), se = F) +
#   facet_wrap(~year(date))
# 
# load("results/mif_tests_sir_sde_CEESBITTBANE009.rda")
# 
# liks <- read_csv("results/loglik_exploration.csv", col_names = F) %>% 
#   set_colnames(c("loglik", names(coef(mf)), "department"))
# 
# 
# # scraps ------------------------------------------------------------------
# 
# 
# logliks <- read_csv("results/loglik_exploration_2.csv", col_names = F)
# load("results/mif_tests_sir_sde_CEESBITTBANE009.rda")
# colnames(logliks) <- c("loglik", "loglik.se", names(coef(mfl[[1]])), "department")
# 
# 
# logliks %>% 
#   filter(loglik > -300) %>% 
#   keep(~sd(.) > 1e-4 | is.character(.)) %>% 
#   group_by(department) %>% 
#   mutate(loglik_std = min(abs(loglik), na.rm = T)/abs(loglik)) %>% 
#   gather(param, val, -contains("loglik"), -department) %>% 
#   ggplot(aes(x = val, color = department)) + 
#   geom_density(aes(y=..scaled..)) +
#   guides(color = "none") +
#   facet_wrap(~param, scales = "free")
# 
# library(GGally)
# p <- ggpairs(logliks %>% 
#                select(loglik, one_of(names(coef(mfl[[1]])))) %>% 
#                keep(~sd(.) > 1e-4) %>% 
#                map_df(~ map_dbl(., ~ifelse(. == 0, NA, .))) %>%
#                bind_cols(logliks[,"department"]),
#              aes(color = department, alpha = I(0.4)),
#              upper = list(continuous = "points", combo = "box_no_facet", discrete = "facetbar", na = "na"),
#              lower = list(continuous = "points", combo = "facethist", discrete = "facetbar", na = "na") ,
#              legend = NULL
# ) +
#   theme_few() +
#   theme(
#     panel.spacing = unit(.5, "lines"),
#     panel.border = element_rect(color = "black", fill = NA, size = .3)
#   ) 

# filtered %>%
#   filter(var %in% c("S", "I", "R")) %>% 
#   select(time, var, q50) %>% 
#   spread(var, q50) %>% mutate(H = I+R+S) %>% 
#   ggplot(aes(x = time, y = H)) + geom_line()
