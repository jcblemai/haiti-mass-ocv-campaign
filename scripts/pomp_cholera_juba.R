# Title: POMP modelling of Cholear epidemic in Juba
# Description: Attempte to fit the model through multiple particle filtering
# Date: Thu Jul  5 09:12:52 2018
# Author: javier.perezsaeez@epfl.ch

# Warning: dt is fraction of a year ! -> so put rates accordily

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

departement <- 'Artibonite'

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

# Load data ---------------------------------------------------------------

# cholera case data from the 2014-2015 epidemic in Juba (South Soudan)
#cases <- read_csv("haiti-data/fromAzman/cases.csv") %>%  
#select('date', departement) %>% 
#     mutate(date = as.Date(date, format = "%Y-%m-%d"),
#            time = dateToYears(date))
# Javier says the second is better
cases <- read_csv("haiti-data/fromAzman/cases.csv")  %>% 
  gather(dep, cases, -date) %>% 
  filter(dep == "Artibonite") %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date))

# get the time of the first datapoint (use %>% filter(time > 2015) to constraint)
t_first_datapnt <- cases %>% slice(1) %>% .[["time"]]

# Estimates of daily rainfall 
#rain <- read_csv("haiti-data/fromAzman/rainfall.csv") %>%  
#select('date', departement) %>%
#mutate(date = as.Date(date, format = "%Y-%m-%d"),
#   time = dateToYears(date)) 

# value of maximal event OK ^
#max_rain <- rain %>%
#filter(year(date)==2015 & month(date) < 10) %>% 
#select(departement) %>% 
#max()


# standardize rainfall  TODO 
#rain %<>% mutate(rain_std = rain/max_rain)


rain <- read_csv("haiti-data/fromAzman/rainfall.csv")  %>% 
  gather(dep, rain, -date) %>% 
  group_by(dep) %>% 
  mutate(max_rain = max(rain), rain_std = rain/max_rain) %>%
  ungroup() %>% 
  filter(dep == "Artibonite") %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date))


make_plots <- F
if(make_plots) {
  # plot the data
  p.data <- ggplot(cases, aes(x = time, y = cases)) + 
    geom_line() + 
    geom_bar(data = rain, aes(y = rain_std  * 100), stat = "identity", fill = "blue") +
    scale_x_continuous(breaks = c(2014, 2015)) 
  
  p.data
}


# Model specification -----------------------------------------------------

## state variable names:
# S:   Susceptibles
# I:   Exposed Infected
# A:   Exposed Asymptomatics (infected but still asymptomatic)
# RI1: Recovered from exposed infected
# RI2: Recovered from exposed infected
# RI3: Recovered from exposed infected
# RA1: Recovered from exposed asymptomatics
# RA2: Recovered from exposed asymptomatics
# RA3: Recovered from exposed asymptomatics
# B:   Bacterial concentration in the environment
# C:   New cases since last observation. Sometime peaks -> data missing and accumulate for two weeks 
#[should add NA to data]. Les cases sont une observation de C par le rmeas.

## data names:
# cases: reported suspected cholera cases (weekly)

## covariate names:
# rain: estimated weekly precipitation standardized with the maximal observation

## parameter names for a departement:
### Pop dynamics
# H:        total population
# D:        density
# gammaI:   recovery from infection
# gammaA:   recovery from asymptomatcs
# mu:       natural mortylity
# alpha:    cholera-induced mortality
# rhoI:      loss of aquired immunity for infected
# rhoA:      loss of aquired immunity for asymptomatics
### Infection dynamics
#### Basic SIRB
# sigma:    sympotmatic to asymptomatic ratio
# betaB:   envionmental (indirect) transmission coefficient
### Bacterial pop dynamics
# mu_B:     bacterial mortality
# thetaA:    bacterial output per asymtomatic person
# thetaI:    bacterial output per infected person
#### Rainfall effects
# lambda: multiplicative rainfall effect on contamination
# r:  power of rainfall effect on contamination
#### Extra-demographic stochasticity
# std_W:    standard deviation of the weiner process to perturb the foi       
### Vaccination
# r_v:      rate of vaccination
# eff_v:    vaccine efficacy
# rho_v:    loss of vaccine-induced immunity
### Measurement model
# epsilon:  under-rerpoting fraction

# Set variables -----------------------------------------------------------

# define stat variable names OK
state_names <- c("S", "I", "A", "RI1", "RI2", "RI3", "RA1", "RA2", "RA3", "B", "C", "W")

# define parameter names for pomp
## process model parameters names to estimate OK
param_proc_est_names <- c("sigma", "betaB", "r", "mu_B", "thetaA", "thetaI", "lambda", "rhoA", "rhoI", "std_W", "epsilon","k")

## initial value parameters to estimate OK
param_iv_est_names <- c("RI1_0")

## fixed process model parameters  OK
param_proc_fixed_names <- c("H", "D", "mu", "alpha", "gammaI", "gammaA")

## fixed initial value parameters OK
param_iv_fixed_names <- c("I_0","A_0", "B_0", "RI2_0", "RI3_0", 
                          "RA1_0", "RA2_0", "RA3_0")

# all paramter names to estimate OK
param_est_names <- c(param_proc_est_names, param_iv_est_names)
# all fixed parameters OK
param_fixed_names <- c(param_proc_fixed_names, param_iv_fixed_names)

# all param names OK
param_names <- c(param_est_names, param_fixed_names)

# names of parameters that are rates MAYBE
param_rates_in_days_names <- c("mu", "alpha", "gammaI", "gammaA", "rhoI", "rhoA")

# names of rate parameters WUT
param_rate_names <- param_names[!str_detect(param_names, "_0|H|k|epsilon|eff|t_|std_W|sigma|alpha_|lambda")] #???


# declare matrix in C for the recoveries in 2014
#cases_2014 <- cases %>% filter(year(date) == 2014)
#cases2014.string <- foreach(r = iter(cases_2014, by = "row"),
#    .combine = c) %do% {
#  sprintf(" {%f, %f} ", r$time, r$cases)
#  } %>% 
#  str_c(collapse = ", \n")

#  matrix_cases2014.string <- str_c(sprintf("double cases2014[%i][%i] = {\n", nrow(cases_2014), 2),
#   cases2014.string,
#   " \n };")


# Measurement model  -------------------------------------------------------

## density

## NegBinomial density (if k -> inf then becomes Poisson) OK TODO S
dmeas <- Csnippet("
  double mean_cases = epsilon * C;
  if (ISNA(cases)) {
    lik = (give_log) ? 0 : 1;
    } else {
        if (S < 10000) {
          lik = (give_log) ? -99999 : 1.0e-18;
          } else {
              lik = dnbinom_mu(cases, k, mean_cases, give_log) ;
          }
      }
      ")

## NegBinomial simulator OK
rmeas <- Csnippet("
  double mean_cases = epsilon * C;
  // cases = mean_cases;
  cases = rnbinom_mu(k, mean_cases);
  ")

# Process model ----------------------------------------------------------------- OK

sirb.rproc <- Csnippet("
  double foi, foi_stoc; // force of infection and its stochastic version
  double dw;            // extra-demographic stochasticity on foi
  double dB;            // deterministic forward time difference of bacteria in the environment
  double k1, k2, k3, k4;  // coefficients of  the Runge-Kutta method
  double r_v_wdn;       // rate of vaccination: 0 if out of time window, r_v if not
  double rate[19];      // vector of all rates in model
  double dN[19];        // vector of transitions between classes during integration timestep

  // force of infection
  foi = betaB * (B / (1 + B));

  if(std_W > 0.0) {
    // white noise (extra-demographic stochasticity)
    dw = rgammawn(std_W, dt);
    // apply stochasticity
    foi_stoc = foi * dw/dt;
    } else {
        foi_stoc = foi;
    }

    // vaccination window
    r_v_wdn = 0.0;

    // define transition rates for each type of event
    // S compartment
    rate[0] = sigma * foi_stoc;   // infections
    rate[1] = (1 - sigma) * foi_stoc;   // asymptomatic infections
    // I compartment
    rate[2] = mu;         // natural deaths
    rate[3] = alpha;      // cholera-induced deaths
    rate[4] = gammaI;      // recovery from infection
    // A compartment (not in order because was added after initial model formulation)
    rate[5] = mu;        // natural death
    rate[6] = gammaA;       // symptoms development
    // RI1,2,3 compartment
    rate[7] = 3*rhoI;        // loss of natural immunity
    rate[8] = mu;         // natural death
    // RI2 compartment
    rate[9] = 3*rhoI;        // loss of natural immunity
    rate[10] = mu;
    // RI3 compartment
    rate[11] = 3*rhoI;        // loss of natural immunity
    rate[12] = mu;
    // RA1,2,3 compartment
    rate[13] = 3*rhoA;        // loss of natural immunity
    rate[14] = mu;         // natural death
    // RA2 compartment
    rate[15] = 3*rhoA;        // loss of natural immunity
    rate[16] = mu;
    // RA3 compartment
    rate[17] = 3*rhoA;        // loss of natural immunity
    rate[18] = mu;


    // simulate all transitions
    reulermultinom(2, S, &rate[0], dt, &dN[0]);
    reulermultinom(3, I, &rate[2], dt, &dN[2]);
    reulermultinom(2, A, &rate[5], dt, &dN[5]);
    reulermultinom(2, RI1, &rate[7], dt, &dN[7]);
    reulermultinom(2, RI2, &rate[9], dt, &dN[9]);
    reulermultinom(2, RI3, &rate[11], dt, &dN[11]);
    reulermultinom(2, RA1, &rate[13], dt, &dN[13]);
    reulermultinom(2, RA2, &rate[15], dt, &dN[15]);
    reulermultinom(2, RA3, &rate[17], dt, &dN[17]);


    // bacteria as continous state variable
    // implement Runge-Kutta integration assuming S, I, R, V* stay constant during dt
    k1 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambda, rain, r, D);
    k2 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambda, rain, r, D);
    k3 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambda, rain, r, D);
    k4 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambda, rain, r, D);
    // bacteria increment
    dB = (k1 + 2*k2 + 2*k3 + k4) / 6.0;

    // update state variables

    I   += dN[0] - dN[2] - dN[3] - dN[4];
    A   += dN[1] - dN[5] - dN[6];
    RI1 += dN[4] - dN[7] - dN[8];
    RI2 += dN[7] - dN[9] - dN[10];
    RI3 += dN[9] - dN[11] - dN[12];
    RA1 += dN[6] - dN[13] - dN[14];
    RA2 += dN[13] - dN[15] - dN[16];
    RA3 += dN[15] - dN[17] - dN[18];
    C   +=  dN[0];
    W   +=  (dw - dt)/std_W;  // standardized i.i.d. white noise
    B += (((dB) < -B) ? (-B + 1.0e-3) : (dB)); // condition to ensure B>0

    // susceptibles so as to match total population
    S = nearbyint(H - I - A - RI1 - RI2 - RI3 - RA1 - RA2 - RA3);
    ")

# C function to compute the time-derivative of bacterial concentration OK
derivativeBacteria.c <- " double fB(int I, int A, double B, 
    double mu_B, double thetaI, double thetaA, double lambda, double rain, double r, double D) {
  double dB;
  dB = -mu_B * B +  (1 + lambda * pow(rain, r)) * D * (thetaI * (double) I + thetaA * (double) A);
  return(dB);
};
"

# C function to compute the initial number of recovered at the start of the simulations given sigma and epsilon.
# computeRecovered2014.c <- "int computeRecovered(double t0, double R_0_2014,  int n_cases2014, double cases2014[][2], double sigma, double rho, double epsilon){
#   double R_0_2015 = 0;
#   // loop over reported cases in 2014 and compute remaning in 205
#   for(int i = 0; i < n_cases2014; i++){
#     R_0_2015 += cases2014[i][1] * (1-sigma)/sigma/epsilon  * exp((cases2014[i][0] - t0) * rho);
# }
# // add the calibrated IC for R in the beginning of 2014 
# R_0_2015 += R_0_2014 * exp((cases2014[0][0] - t0) * rho);

# return(nearbyint(R_0_2015));
# };
# "


# Initializer -------------------------------------------------------------

# initalizeStates <- Csnippet("
#   double m = H ;// /(S_0+E_0+I_0+R_0);
#   A   = nearbyint(A_0 * m/epsilon);
#   I   = nearbyint(I_0 * m/epsilon);
#   R   = computeRecovered(t0, R_0 * H, n_cases2014, cases2014, sigma, rho, epsilon);
#   R   = ((R >= H) ? (m - E - I - 100.0) : (R)); // remove 100 so that S > 0 if the predicted R > H
#   S   = nearbyint(m - E - I - R);
#   B   = 2.0/epsilon * theta/mu_B; // custom initial conditions equivalent to the 'forcing' in the continous model
#   C   = 0;
#   W   = 0;
#   Vtot  = 0;
#   ")

# TODO
initalizeStates <- Csnippet("
  double m = H ;// /(S_0+E_0+I_0+R_0);
  A   = 0;
  I   = 0;
  RI1   = 100;
  RI2   = 100;
  RI3   = 100;
  RA1   = 100;
  RA2   = 100;
  RA3   = 100;
  S   = nearbyint(m - A - I - RI1 - RI2 - RI3 - RA1 -RA2 -RA3);
  B   = 2.0/epsilon * thetaI/mu_B; // TODO custom initial conditions equivalent to the 'forcing' in the continous model
  C   = 0;
  W   = 0;
  ")


# Parameter transformations -----------------------------------------------
# use log for positive parameters and logit for parmaters in [0,1]

toEstimationScale <- Csnippet("
  Tsigma = logit(sigma);
  TbetaB = log(betaB);
  Tmu_B = log(mu_B);
  TthetaA = log(thetaA);
  TthetaI = log(thetaI);
  TrhoI = log(rhoI);
  TrhoA = log(rhoA);
  Tlambda = log(lambda);
  Tr = log(r);
  Tstd_W = log(std_W);
  Tepsilon = log(epsilon);
  Tk = log(k);
  TRI1_0 = logit(RI1_0);
  TB_0 = log(B_0);
  ")

fromEstimationScale <- Csnippet("
  Tsigma = expit(sigma);
  TbetaB = exp(betaB);
  Tmu_B = exp(mu_B);
  TthetaI = exp(thetaI);
  TthetaA = exp(thetaA);
  TrhoA = exp(rhoA);
  TrhoI = exp(rhoI);
  Tlambda = exp(lambda);
  Tr = exp(r);
  Tstd_W = exp(std_W);
  Tepsilon = exp(epsilon);
  Tk = exp(k);
  TRI1_0 = expit(RI1_0);
  TB_0 = exp(B_0);
  ")

# Build pomp object -------------------------------------------------------

# input parameters to the model
input_parameters <- yaml::read_yaml("data/input_parameters.yaml")

# Start and end dates of epidemic
t_start <- dateToYears(as.Date(input_parameters$t_start))
t_end <- dateToYears(as.Date(input_parameters$t_end))

# get fixed process paramteres to input
fixed_input_parameters <- as_vector(input_parameters[map_lgl(names(input_parameters), ~ . %in% param_fixed_names)])

# set fixed process parameters
param_proc_fixed <- set_names(seq_along(param_proc_fixed_names) * 0, param_proc_fixed_names)
param_proc_fixed[names(fixed_input_parameters)] <- fixed_input_parameters

# Initialize the fixed parameters
param_fixed <-  set_names(seq_along(param_fixed_names) * 0, param_fixed_names)
param_fixed[param_proc_fixed_names] <- as.numeric(param_proc_fixed)

# Initial Conditions based on forcing WUT
param_fixed["A_0"] <- 3 / param_fixed["H"]
param_fixed["I_0"] <- 2 / param_fixed["H"]
param_fixed["B_0"] <- 0 # B0 depends on epsilon and sigma

# Initialize the parameters to estimate (just initial guesses)
param_est <- set_names(seq_along(param_est_names) * 0, param_est_names)
param_est["sigma"] <- .2
param_est["rhoA"] <- 1/(365*3)
param_est["rhoI"] <- 1/(365*3)
param_est["betaB"] <- .1
param_est["mu_B"] <-  365/5
param_est["thetaA"] <- .01
param_est["thetaI"] <- .01
param_est["lambda"] <- 100
param_est["r"] <- 3
param_est["std_W"] <- .001
param_est["epsilon"] <- .5
param_est["k"] <- 10001
param_est["RI1_0"] <- 0.1

# rate of simulation in fractions of years
dt_yrs <- 1/365.25 * .1

# adjust the rate parameters depending on the integration delta time in years (some parameter inputs given in days) TODO CHECK
params <- c(param_est, param_fixed)
params[param_rates_in_days_names] <- params[param_rates_in_days_names] * 365.25

sirb_cholera <- pomp(
  # set data
  data = cases %>% 
    filter(time > t_start & time < (t_end + 0.01)) %>% 
    select(time, cases) , 
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
  # covariates
  covar = rain  %>% 
    filter(time > (t_start - 0.01) & time < (t_end + 0.01)) %>% 
    select(time, rain_std) %>% 
    rename(rain = rain_std),
  tcovar = "time",
  # names of state variables
  statenames = state_names,
  # names of accumulator variables to be re-initalized at each observation timestep 
  # (C for cases, W for the white noise just for plotting)
  zeronames = c("C", "W"),
  # names of paramters
  paramnames = param_names,
  # names of covariates
  covarnames = "rain",
  initializer = initalizeStates,
  toEstimationScale = toEstimationScale,
  fromEstimationScale = fromEstimationScale,
  # global C definitions
  globals = str_c(
    #    sprintf("double t_vacc_start = %f; double t_vacc_end = %f;", t_vacc_start, t_vacc_end),
    #    sprintf("double t0 = %f;",  t_first_datapnt - dt_yrs),
    derivativeBacteria.c,
    #matrix_cases2014.string, 
    #sprintf("int n_cases2014 = %i;",  nrow(cases_2014)),
    #computeRecovered2014.c,
    sep = " ")
)

# save pomp object for further use
save(sirb_cholera, file = "data/sirb_cholera_pomped.rda")


p <- simulate(sirb_cholera, nsim = 10, as.data.frame = T) %>% 
  gather(variable, value, -time, -rain, -sim)  %>% 
  ggplot(aes(x = time, y = value, color = sim)) + 
  geom_line() + 
  facet_wrap(~variable, scales = "free_y")
# use coef(sirb_cholera)["mu_B"] <- 365/5 to change a parameter

print(p)