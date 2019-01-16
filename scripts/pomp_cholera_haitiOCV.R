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
output_dir <- "output/"


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  # default departement
  args[1] = "Sud"
}
departement <- args[1]

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

# input parameters to the model
input_parameters <- yaml::read_yaml("haiti-data/input_parameters.yaml")

# Start and end dates of epidemic
t_start <- dateToYears(as.Date(input_parameters$t_start))
t_end <- dateToYears(as.Date(input_parameters$t_end))


cases <- read_csv("haiti-data/fromAzman/cases_corrected.csv")  %>% 
  gather(dep, cases, -date) %>% 
  filter(dep == departement) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date))

case_dates <- with(cases %>% 
                       filter(time > t_start - 0.01 & time < (t_end + 0.01)),
                   seq.Date(min(date), max(date), by = "1 week")
                   )
# check if weekly reports are missing
missing_dates <- setdiff(case_dates, cases$date) %>% as.Date(origin = as.Date("1970-01-01"))

if (length(missing_dates) > 0) {
# fill in the data
  cases %<>% 
    bind_rows(slice(cases, 1:2) %>% 
                mutate(date = missing_dates,
                      cases = NA,
                      time = dateToYears(date))) %>%
    arrange(time)
}

rain <- read_csv("haiti-data/fromAzman/rainfall.csv")  %>% 
  gather(dep, rain, -date) %>% 
  group_by(dep) %>% 
  ungroup() %>% 
  filter(dep == departement) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date)) %>%
  filter(time > t_start - 0.01 & time < (t_end + 0.01)) %>%
  mutate(max_rain = max(rain), rain_std = rain/max_rain) 


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
state_names <- c("S", "I", "A", "RI1", "RI2", "RI3", "RA1", "RA2", "RA3",
 "VSd", "VRI1d", "VRI2d", "VRI3d", "VRA1d", "VRA2d", "VRA3d",
 "VSdd", "VRI1dd", "VRI2dd", "VRI3dd", "VRA1dd", "VRA2dd", "VRA3dd",
"VSd_alt", "VRI1d_alt", "VRI2d_alt", "VRI3d_alt", "VRA1d_alt", "VRA2d_alt", "VRA3d_alt",
 "VSdd_alt", "VRI1dd_alt", "VRI2dd_alt", "VRI3dd_alt", "VRA1dd_alt", "VRA2dd_alt", "VRA3dd_alt",
"B", "C", "W")



# define parameter names for pomp
## process model parameters names to estimate OK
param_proc_est_names <- c("sigma", "betaB", "mu_B", "thetaI", "XthetaA",  "lambdaR", "r",  
                          "gammaI", "gammaA", "rhoA", "XrhoI", "foi_add", "epsilon","k","std_W")

## initial value parameters to estimate OK
param_iv_est_names <- c("Rtot_0")

## fixed process model parameters  OK
param_proc_fixed_names <- c("H", "D", "mu", "alpha")

# Vaccination senario:
param_vacc_fixed_names <- c("t_vacc_start", "t_vacc_end", "p1d_reg", "r_v_year")

# all paramter names to estimate OK
param_est_names <- c(param_proc_est_names, param_iv_est_names)
# all fixed parameters OK
param_fixed_names <- c(param_proc_fixed_names, param_vacc_fixed_names)
# all param names OK
param_names <- c(param_est_names, param_fixed_names)

# names of parameters that are rates (MAYBE) (because time 365 since timestep is year) r_v shoudl be here
param_rates_in_days_names <- c("mu", "alpha", "gammaI", "gammaA", "rhoA") #muB

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

sirb_file <- 'scripts/sirb_model_vacc.c'

sirb.rproc <- Csnippet(readChar(sirb_file, file.info(sirb_file)$size))

# C function to compute the time-derivative of bacterial concentration OK
derivativeBacteria.c <- " double fB(int I, int A, double B, 
    double mu_B, double thetaI, double XthetaA, double lambdaR, double rain, double r, double D) {
  double thetaA = thetaI * XthetaA;
  double dB;
  dB = -mu_B * B +  (1 + lambdaR * pow(rain, r)) * D * (thetaI * (double) I + thetaA * (double) A);
  return(dB);
};
"


# Initializer -------------------------------------------------------------
initalizeStates <- Csnippet("
  A     = nearbyint((1-sigma)/sigma  * 1/epsilon * cases0/7 * 365 /(mu+gammaA));
  I     = nearbyint(1/epsilon * cases0/7 * 365 /(mu+alpha+gammaI))  ;  // Steady state
  RI1   = nearbyint(sigma * Rtot_0*H/3.0);
  RI2   = nearbyint(sigma * Rtot_0*H/3.0);
  RI3   = nearbyint(sigma * Rtot_0*H/3.0);
  RA1   = nearbyint((1-sigma) * Rtot_0*H/3.0);
  RA2   = nearbyint((1-sigma) * Rtot_0*H/3.0);
  RA3   = nearbyint((1-sigma) * Rtot_0*H/3.0);
  if (A + I + RI1 + RI2 + RI3 + RA1 + RA2 + RA3 >= H)
  {
    double R_tot = H - A - I - 100.0;
    if (R_tot <= 0)
    {
      I     = nearbyint(H - 100);
      A     = nearbyint(0);
      R_tot = nearbyint(0);
    }
    RI1   = nearbyint(sigma * R_tot/3.0);
    RI2   = nearbyint(sigma * R_tot/3.0);
    RI3   = nearbyint(sigma * R_tot/3.0);
    RA1   = nearbyint((1-sigma) * R_tot/3.0);
    RA2   = nearbyint((1-sigma) * R_tot/3.0);
    RA3   = nearbyint((1-sigma) * R_tot/3.0);
  }
  S   = nearbyint(H - A - I - RI1 - RI2 - RI3 - RA1 - RA2 - RA3);
  B   = 2.0/epsilon * thetaI/mu_B; // TODO custom initial conditions equivalent to the 'forcing' in the continous model
  C   = 0;
  W   = 0;
  VSd = 0;
  VRI1d = 0;
  VRI2d = 0;
  VRI3d = 0;
  VRA1d = 0;
  VRA2d = 0;
  VRA3d = 0;
  VSdd = 0;
  VRI1dd = 0;
  VRI2dd = 0;
  VRI3dd = 0;
  VRA1dd = 0;
  VRA2dd = 0;
  VRA3dd = 0;
  VSd_alt = 0;
  VRI1d_alt = 0;
  VRI2d_alt = 0;
  VRI3d_alt = 0;
  VRA1d_alt = 0;
  VRA2d_alt = 0;
  VRA3d_alt = 0;
  VSdd_alt = 0;
  VRI1dd_alt = 0;
  VRI2dd_alt = 0;
  VRI3dd_alt = 0;
  VRA1dd_alt = 0;
  VRA2dd_alt = 0;
  VRA3dd_alt = 0;
   ")


eff_v.c <- paste0(readChar('scripts/v_eff.c', file.info(sirb_file)$size), " double eff_v_1d(double t_since_vacc, int scenario) {
  if (t_since_vacc < 1) 
    return eff_v_2d(t_since_vacc, scenario);
  else
   return 0;
};
")




# Parameter transformations -----------------------------------------------
# use log for positive parameters and logit for parmaters in [0,1]

toEstimationScale <- Csnippet("
  Tsigma = logit(sigma);
  TbetaB = log(betaB);
  Tmu_B = log(mu_B);
  TXthetaA = logit(XthetaA);
  TthetaI = log(thetaI);
  TXrhoI = logit(XrhoI);
  TrhoA = log(rhoA);
  TlambdaR = log(lambdaR);
  Tr = log(r);
  Tstd_W = log(std_W);
  Tepsilon = logit(epsilon);
  Tk = log(k);
  TRtot_0 = logit(Rtot_0);
  Tfoi_add = log(foi_add);
  TgammaA = log(gammaA);
  TgammaI = log(gammaI);
  ")

fromEstimationScale <- Csnippet("
  Tsigma = expit(sigma);
  TbetaB = exp(betaB);
  Tmu_B = exp(mu_B);
  TthetaI = exp(thetaI);
  TXthetaA = expit(XthetaA);
  TrhoA = exp(rhoA);
  TXrhoI = expit(XrhoI);
  TlambdaR = exp(lambdaR);
  Tr = exp(r);
  Tstd_W = exp(std_W);
  Tepsilon = expit(epsilon);
  Tk = exp(k);
  TRtot_0 = expit(Rtot_0);
  Tfoi_add = exp(foi_add);
  TgammaA = exp(gammaA);
  TgammaI = exp(gammaI);
  ")

# Build pomp object -------------------------------------------------------

# input parameters to the model
input_parameters <- yaml::read_yaml("haiti-data/input_parameters.yaml")

# Start and end dates of epidemic
t_start <- dateToYears(as.Date(input_parameters$t_start))
t_end <- dateToYears(as.Date(input_parameters$t_end))

# get fixed process paramteres to input
fixed_input_parameters <- as_vector(input_parameters[map_lgl(names(input_parameters), ~ . %in% param_fixed_names)])


# set fixed process parameters
param_proc_fixed <- set_names(seq_along(param_proc_fixed_names) * 0, param_proc_fixed_names)
param_proc_fixed[names(fixed_input_parameters)] <- fixed_input_parameters

populations  <- unlist(flatten(input_parameters["population"]))
densities <- unlist(flatten(input_parameters["density"]))
param_proc_fixed['H'] <- populations[departement]
param_proc_fixed['D'] <- densities[departement]

p1d_alt_year  <- unlist(flatten(input_parameters["p1d_alt_year"]))
nb_doses_alt_year <- unlist(flatten(input_parameters["nb_doses_alt_year"]))
t_vacc_start_alt  <- unlist(flatten(input_parameters["t_vacc_start_alt"]))
t_vacc_end_alt <- unlist(flatten(input_parameters["t_vacc_end_alt"]))

# Vaccination information:
t_vacc_start_alt = dateToYears(as.Date(t_vacc_start_alt[departement]))
t_vacc_end_alt   = dateToYears(as.Date(t_vacc_end_alt[departement]))
r_v_alt_year = nb_doses_alt_year[departement]/(t_vacc_end_alt - t_vacc_start_alt)
p1d_alt = p1d_alt_year[departement]



# Initialize the fixed parameters
param_fixed <-  set_names(seq_along(param_fixed_names) * 0, param_fixed_names)
param_fixed[param_proc_fixed_names] <- as.numeric(param_proc_fixed)  # Does not work for gammaI TODO


#Cases in the last report:
# declare matrix in C for the infected before the strat date in 2014
cases_at_t_start <- cases %>% filter(dateToYears(date) <= t_start) %>% tail(n=1)%>% select('cases') %>% unlist()
cases_at_t_start.string <- sprintf("double cases0 = %i;", cases_at_t_start)


# Initialize the parameters to estimate (just initial guesses)
param_est <- set_names(seq_along(param_est_names) * 0, param_est_names)
param_est["sigma"] <- .2
param_est["rhoA"] <- 1/(365*3)
param_est["XrhoI"] <- 0.5
param_est["betaB"] <- .001
param_est["mu_B"] <-  365/5
param_est["XthetaA"] <- 0.5
param_est["thetaI"] <- .01
param_est["lambdaR"] <- 10
param_est["r"] <- 1
param_est["std_W"] <- .001
param_est["epsilon"] <- .5
param_est["k"] <- 10001
param_est["Rtot_0"] <- 0.35
param_est["foi_add"] <- 0.0
param_est["gammaA"] <- 1/5
param_est["gammaI"] <- 1/5

# rate of simulation in fractions of years
dt_yrs <- 1/365.25 * .2

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
    sprintf("double t_vacc_start_alt = %f; double t_vacc_end_alt = %f;", t_vacc_start_alt, t_vacc_end_alt),
    sprintf("double r_v_alt = %f;", r_v_alt_year),
    sprintf("double p1d_alt = %f;", p1d_alt),
    derivativeBacteria.c,
    cases_at_t_start.string,
    eff_v.c,
    sep = " ")
)

# save pomp object for further use
save(sirb_cholera, file = paste0(output_dir, departement, "/sirb_cholera_pomped_", departement, ".rda"))
