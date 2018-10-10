# Title: POMP modelling of Cholear epidemic in Juba
# Description: Attempte to fit the model through multiple particle filtering
# Date: Thu Jul  5 09:12:52 2018
# Author: javier.perezsaeez@epfl.ch



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

# Load data ---------------------------------------------------------------

# cholera case data from the 2014-2015 epidemic in Juba (South Soudan)
cases <- read_csv("data/case_data.csv") %>% 
  mutate(date = as.Date(date, format = "%d-%b-%y"),
         time = dateToYears(date))

# get the time of the first datapoint in 2015 (used for simulations)
t_first_datapnt <- cases %>% filter(time > 2015) %>% slice(1) %>% .[["time"]]

# Estimates of daily rainfall 
rain <- read_csv("data/rainfall_data.csv") %>% 
  mutate(date = as.Date(date, format = "%d-%b-%y"),
         time = dateToYears(date)) 

# value of maximal event during the 2015 epidemic
max_rain2015 <- rain %>%
  filter(year(date)==2015 & month(date) < 10) %>% 
  select(rain) %>% 
  max()

# standardize rainfall
rain %<>% mutate(rain_std = rain/max_rain2015)

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
# S:  Susceptibles
# I:  Exposed Infected
# A:  Exposed Asymptomatics (infected but still asymptomatic)
# RI1,RI2,RI3:  Recovered from exposed infected
# RA1,RA2,RA3:  Recovered from exposed Asymptomatics
# B:  Bacterial concentration in the environment
# C: Cumulative cases

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
# beta:   envionmental (indirect) transmission coefficient
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

# define stat variable names
state_names <- c("S", "I", "A", "RI1", "RI2", "RI3", "RA1", "RA2", "RA3", "B", "C", "W")#,  "VS", "VE", "VI", "VR", "Vtot")

# define parameter names for pomp
## process model parameters names to estimate
param_proc_est_names <- c("sigma", "beta", "r", "mu_B", "thetaA", "thetaI", "lambda", "gammaI", "gammaA", "rhoA", "rhoI", "std_W", "epsilon","k")

## initial value parameters to estimate
param_iv_est_names <- c("RI1_0","RI2_0","RI3_0","RA1_0","RA2_0","RA3_0")

## fixed process model parameters 
param_proc_fixed_names <- c("H", "D", "mu", "phi", "gamma", "alpha")

## fixed initial value parameters
param_iv_fixed_names <- c("E_0", "I_0","A_0", "B_0", "VS_0", "VE_0", "VI_0", "VR_0")

# all paramter names to estimate
param_est_names <- c(param_proc_est_names, param_iv_est_names)
# all fixed parameters
param_fixed_names <- c(param_proc_fixed_names, param_iv_fixed_names)

# all param names
param_names <- c(param_est_names, param_fixed_names)

# names of parameters that are rates
param_rates_in_days_names <- c("mu", "alpha", "phi", "gamma", "rho", "r_v")

# names of rate parameters
param_rate_names <- param_names[!str_detect(param_names, "_0|H|k|epsilon|eff|t_|std_W|sigma|alpha_|lambda")]


# declare matrix in C for the recoveries in 2014
cases_2014 <- cases %>% filter(year(date) == 2014)
cases2014.string <- foreach(r = iter(cases_2014, by = "row"),
        .combine = c) %do% {
          sprintf(" {%f, %f} ", r$time, r$cases)
        } %>% 
  str_c(collapse = ", \n")

matrix_cases2014.string <- str_c(sprintf("double cases2014[%i][%i] = {\n", nrow(cases_2014), 2),
                                 cases2014.string,
                                 " \n };")


# Measurment model  -------------------------------------------------------

# measurement model
## density

## NegBinomial density (if k -> inf then becomes Poisson)
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

## NegBinomial simulator
rmeas <- Csnippet("
  double mean_cases = epsilon * C;
  cases = rnbinom_mu(k, mean_cases);
")

# Process model -----------------------------------------------------------------

sirb.rproc <- Csnippet("
  double foi, foi_stoc; // force of infection and its stochastic version
  double dw;            // extra-demographic stochasticity on foi
  double dB;            // deterministic forward time difference of bacteria in the environment
  double k1, k2, k3, k4;  // coefficients of  the Runge-Kutta method
  double r_v_wdn;       // rate of vaccination: 0 if out of time window, r_v if not
  double rate[23];      // vector of all rates in model
  double dN[23];        // vector of transitions between classes during integration timestep

  // force of infection
  foi = beta_B * (B / (1 + B)) * (1 + lambda_E * pow(rain, alpha_E)) + beta_I * (I + VI) / H;

  if(std_W > 0.0) {
    // white noise (extra-demographic stochasticity)
    dw = rgammawn(std_W, dt);
    // apply stochasticity
    foi_stoc = foi * dw/dt;
  } else {
    foi_stoc = foi;
  }
  
  // vaccination window
  if (t >= t_vacc_start && t <= (t_vacc_end + dt)) 
    r_v_wdn = (r_v / (S + E + R));
  else 
    r_v_wdn = 0.0;
  
  // define transition rates for each type of event
  // S compartment
  rate[0] = sigma * foi_stoc;   // infections
  rate[1] = (1 - sigma) * foi_stoc;   // asymptomatic infections
  rate[2] = r_v_wdn;    // vaccinations
  // I compartment
  rate[3] = mu;         // natural deaths
  rate[4] = alpha;      // cholera-induced deaths
  rate[5] = gamma;      // recovery from infection
  // R compartment
  rate[6] = rho;        // loss of natural immunity
  rate[7] = mu;         // natural death
  rate[8] = r_v_wdn;    // vaccinations
  // VS compartment
  rate[9] = mu;          // natural death
  rate[10] = sigma * (1 - eff_v) * foi_stoc; // symptomatic infections
  rate[11] = (1 - sigma) * (1 - eff_v) * foi_stoc; // asymptomatic infections
  // VI compartment
  rate[12] = mu;          // natural death
  rate[13] = alpha;       // cholera-induced death
  rate[14] = gamma;       // recovery
  // VR compartment
  rate[15] = mu;          // natural death
  rate[16] = rho_v;       // loss of vaccine immunity

  // E compartment (not in order because was added after initial model formulation)
  rate[17] = mu;        // natural death
  rate[18] = r_v_wdn;   // vaccination
  rate[19] = phi;       // symptoms development
  // VE compartment
  rate[20] = mu;          // natural death
  rate[21] = alpha;       // cholera-induced death
  rate[22] = phi;         // symptoms development

  // simulate all transitions
  reulermultinom(3, S, &rate[0], dt, &dN[0]);
  reulermultinom(3, I, &rate[3], dt, &dN[3]);
  reulermultinom(3, R, &rate[6], dt, &dN[6]);
  reulermultinom(3, VS, &rate[9], dt, &dN[9]);
  reulermultinom(3, VI, &rate[12], dt, &dN[12]);
  reulermultinom(2, VR, &rate[15], dt, &dN[15]);
  reulermultinom(3, E, &rate[17], dt, &dN[17]);
  reulermultinom(3, VE, &rate[20], dt, &dN[20]);

  // bacteria as continous state variable
  // implement Runge-Kutta integration assuming S, I, R, V* stay constant during dt
  k1 = dt * fB(I, VI, B, mu_B, theta, lambda_R, rain, alpha_R);
  k2 = dt * fB(I, VI, B + k1/2, mu_B, theta, lambda_R, rain, alpha_R);
  k3 = dt * fB(I, VI, B + k2/2, mu_B, theta, lambda_R, rain, alpha_R);
  k4 = dt * fB(I, VI, B + k3, mu_B, theta, lambda_R, rain, alpha_R);
  // bacteria increment
  dB = (k1 + 2*k2 + 2*k3 + k4) / 6.0;

  // update state variables
  E   += dN[0] - dN[17] - dN[18] - dN[19];
  I   += dN[19] - dN[3] - dN[4] - dN[5];
  R   += dN[1] - dN[6] - dN[7] - dN[8] + dN[5];
  VS  += dN[2] - dN[9] + dN[16] - dN[10] - dN[11];
  VE  += dN[18] + dN[10] - dN[20] - dN[21] - dN[22];
  VI  += dN[22] - dN[13] - dN[14];
  VR  += dN[8] - dN[15] - dN[16] + dN[14] + dN[11];
  C   +=  dN[19] + dN[22];
  W   +=  (dw - dt)/std_W;  // standardized i.i.d. white noise
  Vtot += dN[2] + dN[8];
  B += (((dB) < -B) ? (-B + 1.0e-3) : (dB)); // condition to ensure B>0


  // susceptibles so as to match total population
  S = nearbyint(H - I - E - R - VS - VE - VI - VR);
")

# C function to compute the time-derivative of bacterial concentration
derivativeBacteria.c <- " double fB(int I, int VI, double B, 
double mu_B, double theta, double lambda_R, double rain, double alpha_R) {
  double dB;
  dB = -mu_B * B + theta * (1 + lambda_R * pow(rain, alpha_R)) * ((double) I + (double) VI);
  return(dB);
};
"

# C function to compute the initial number of recovered at the start of the simulations given sigma and epsilon
computeRecovered2014.c <- "int computeRecovered(double t0, double R_0_2014,  int n_cases2014, double cases2014[][2], double sigma, double rho, double epsilon){
  double R_0_2015 = 0;
  // loop over reported cases in 2014 and compute remaning in 205
  for(int i = 0; i < n_cases2014; i++){
    R_0_2015 += cases2014[i][1] * (1-sigma)/sigma/epsilon  * exp((cases2014[i][0] - t0) * rho);
  }
  // add the calibrated IC for R in the beginning of 2014 
  R_0_2015 += R_0_2014 * exp((cases2014[0][0] - t0) * rho);
  
  return(nearbyint(R_0_2015));
};
"
# Deterministic skeleton --------------------------------------------------
sirb.skeleton <- Csnippet("
                       double foi; // force of infection and its stochastic version
                       double r_v_wdn;       // rate of vaccination: 0 if out of time window, r_v if not
                       double rate[23];      // vector of all rates in model
                       
                       // force of infection
                       foi = beta_B * (B / (1 + B)) * (1 + lambda_E * pow(rain, alpha_E)) + beta_I * (I + VI) / H;
                       
                      // vaccination window
                       if (t >= t_vacc_start && t <= (t_vacc_end + 1/365.25)) 
                       r_v_wdn = (r_v / (S + E + R));
                       else 
                       r_v_wdn = 0;
                       
                       // define transition rates for each type of event
                       // S compartment
                       rate[0] = sigma * foi;   // infections
                       rate[1] = (1 - sigma) * foi;   // asymptomatic infections
                       rate[2] = r_v_wdn;    // vaccinations
                       // I compartment
                       rate[3] = mu;         // natural deaths
                       rate[4] = alpha;      // cholera-induced deaths
                       rate[5] = gamma;      // recovery from infection
                       // R compartment
                       rate[6] = rho;        // loss of natural immunity
                       rate[7] = mu;         // natural death
                       rate[8] = r_v_wdn;    // vaccinations
                       // VS compartment
                       rate[9] = mu;          // natural death
                       rate[10] = sigma * (1 - eff_v) * foi; // symptomatic infections
                       rate[11] = (1 - sigma) * (1 - eff_v) * foi; // asymptomatic infections
                       // VI compartment
                       rate[12] = mu;          // natural death
                       rate[13] = alpha;       // cholera-induced death
                       rate[14] = gamma;       // recovery
                       // VR compartment
                       rate[15] = mu;          // natural death
                       rate[16] = rho_v;       // loss of vaccine immunity
                       
                       // E compartment
                       rate[17] = mu;        // natural death
                       rate[18] = r_v_wdn;   // vaccination
                       rate[19] = phi;       // symptoms development
                       // VE compartment
                       rate[20] = mu;          // natural death
                       rate[21] = alpha;       // cholera-induced death
                       rate[22] = phi;         // symptoms development
                       
                       // update state variables
                       DE  = rate[0] * S - (rate[17] + rate[18] + rate[19]) * E;
                       DI  = rate[19] * E - (rate[3] + rate[4] + rate[5]) * I;
                       DR  = rate[1] * S - (rate[6] + rate[7] + rate[8]) * R + rate[5] * I;
                       DVS = rate[2] * S + rate[16] * VR - (rate[9] + rate[10] + rate[11]) * VS;
                       DVE = rate[18] * E + rate[10] * VS - (rate[20] + rate[21] + rate[22]) * VE;
                       DVI = rate[22] * VE - (rate[13] + rate[14]) * VI;
                       DVR = rate[8] * R - (rate[15] + rate[16]) * VR + rate[14] * VI + rate[11] * VS;
                       DC  = rate[19] * E + rate[22] * VE;
                       DW  = foi;  // standardized i.i.d. white noise
                       DVtot = rate[2] * S + rate[18] * E + rate[8] * R;
                       
                       // bacteria as continous state variable
                       DB = -mu_B * B + theta * (1 + lambda_R * pow(rain, alpha_R)) * (I + VI);

                       // susceptibles so as to match total population
                       DS = -(DE + DI + DR + DVS + DVI + DVE + DVR);
                          ")



# Initializer -------------------------------------------------------------

initalizeStates <- Csnippet("
  double m = H ;// /(S_0+E_0+I_0+R_0);
  E   = nearbyint(E_0 * m/epsilon);
  I   = nearbyint(I_0 * m/epsilon);
  R   = computeRecovered(t0, R_0 * H, n_cases2014, cases2014, sigma, rho, epsilon);
  R   = ((R >= H) ? (m - E - I - 100.0) : (R)); // remove 100 so that S > 0 if the predicted R > H
  S   = nearbyint(m - E - I - R);
  VS  = nearbyint(VS_0 * m);
  VE  = nearbyint(VE_0 * m);
  VI  = nearbyint(VI_0 * m);
  VR  = nearbyint(VR_0 * m);
  B   = 2.0/epsilon * theta/mu_B; // custom initial conditions equivalent to the 'forcing' in the continous model
  C   = 0;
  W   = 0;
  Vtot  = 0;
")


# Parameter transformations -----------------------------------------------
# use log for positive parameters and logit for parmaters in [0,1]

toEstimationScale <- Csnippet("
  Tsigma = logit(sigma);
  Tbeta_B = log(beta_B);
  Tbeta_I = log(beta_I);
  Tmu_B = log(mu_B);
  Ttheta = log(theta);
  Trho = log(rho);
  Tlambda_E = log(lambda_E);
  Tlambda_R = log(lambda_R);
  Talpha_E = log(alpha_E);
  Talpha_R = log(alpha_R);
  Tstd_W = log(std_W);
  Tepsilon = log(epsilon);
  Tk = log(k);
  TR_0 = logit(R_0);
  TB_0 = log(B_0);
")

fromEstimationScale <- Csnippet("
  Tsigma = expit(sigma);
  Tbeta_B = exp(beta_B);
  Tbeta_I = exp(beta_I);
  Tmu_B = exp(mu_B);
  Ttheta = exp(theta);
  Trho = exp(rho);
  Tlambda_E = exp(lambda_E);
  Tlambda_R = exp(lambda_R);
  Talpha_E = exp(alpha_E);
  Talpha_R = exp(alpha_R);
  Tstd_W = exp(std_W);
  Tepsilon = exp(epsilon);
  Tk = exp(k);
  TR_0 = expit(R_0);
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

# Start and end dates of vaccination
t_vacc_start <- dateToYears(as.Date(input_parameters$t_vacc_start))
t_vacc_end <- dateToYears(as.Date(input_parameters$t_vacc_end))

# Initialize the fixed parameters
param_fixed <-  set_names(seq_along(param_fixed_names) * 0, param_fixed_names)
param_fixed[param_proc_fixed_names] <- as.numeric(param_proc_fixed)

# Initial Conditions based on forcing
param_fixed["E_0"] <- 3 / param_fixed["H"]
param_fixed["I_0"] <- 2 / param_fixed["H"]
param_fixed["B_0"] <- 0 # B0 depends on epsilon and sigma

# Initialize the parameters to estimate (just initial guesses)
param_est <- set_names(seq_along(param_est_names) * 0, param_est_names)
param_est["sigma"] <- .2
param_est["rho"] <- 1/(365*3)
param_est["beta_B"] <- 3
param_est["beta_I"] <- 1
param_est["mu_B"] <-  4000
param_est["theta"] <- 1
param_est["lambda_E"] <- 0
param_est["lambda_R"] <- 0
param_est["alpha_E"] <- 1
param_est["alpha_R"] <- 1
param_est["std_W"] <- .001
param_est["epsilon"] <- .5
param_est["k"] <- 1
param_est["R_0"] <- 0.1

# rate of simulation in fractions of years
dt_yrs <- 1/365.25 * .1

# adjust the rate parameters depending on the integration delta time in years (some parameter inputs given in days)
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
  # deterministic skeleton
  skeleton = vectorfield(sirb.skeleton),
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
    sprintf("double t_vacc_start = %f; double t_vacc_end = %f;", t_vacc_start, t_vacc_end),
    sprintf("double t0 = %f;",  t_first_datapnt - dt_yrs),
    derivativeBacteria.c,
    matrix_cases2014.string, 
    sprintf("int n_cases2014 = %i;",  nrow(cases_2014)),
    computeRecovered2014.c,
    sep = " ")
)

# save pomp object for further use
save(sirb_cholera, file = "data/sirb_cholera_pomped.rda")
 
