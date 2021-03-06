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
  args[1] = "Artibonite"
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

#if (length(missing_dates) > 0) {
## fill in the data
#  cases %<>% 
#    bind_rows(slice(cases, 1:2) %>% 
#                mutate(date = missing_dates,
#                      cases = NA,
#                      time = dateToYears(date))) %>%
#    arrange(time)
#} # forgot how that work ^ and cannot adapt it to new data

rain <- read_csv("haiti-data/fromAzman/rainfall.csv")  %>% 
  gather(dep, rain, -date) %>% 
  group_by(dep) %>% 
  ungroup() %>% 
  filter(dep == departement) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date)) %>%
  filter(time > t_start - 0.01 & time < (t_end + 0.01)) %>%
  mutate(max_rain = max(rain), rain_std = rain/max_rain)

cases_other_dept <- read_csv("haiti-data/fromAzman/cases_corrected.csv")  %>% 
  gather(dep, cases, -date) %>% 
  filter(dep != departement) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date))

cases_other_dept <- aggregate(cases_other_dept$cases, by=list(Category=cases_other_dept$time), FUN=sum, na.rm=TRUE, na.action=NULL) %>%
  mutate(time = Category) %>%
  mutate(cases = x)

cases_covar <- cases_other_dept
#cases_covar <- read_csv("covar_zero.csv")  %>% 
#  gather(dep, cases, -date) %>% 
#  filter(dep != departement) %>% 
#  mutate(date = as.Date(date, format = "%Y-%m-%d"),
#         time = dateToYears(date))

#cases_covar <- aggregate(cases_covar$cases, by=list(Category=cases_covar$time), FUN=sum, na.rm=TRUE, na.action=NULL) %>%
#  mutate(time = Category) %>%
#  mutate(cases = x)


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
# cas_def: parameter for the change of reporting induced by the CDC cases definition on January 1, 2018
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
                          "gammaI", "gammaA", "rhoA", "XrhoI", "foi_add", "epsilon","k","std_W", "cas_def")

## initial value parameters to estimate OK
param_iv_est_names <- c("Rtot_0")

## fixed process model parameters  OK
param_proc_fixed_names <- c("H", "D", "mu", "alpha", "cases_ext")

# Vaccination senario:
param_vacc_fixed_names <- c("t_vacc_start", "t_vacc_end", "p1d_reg", "r_v_year")

# all paramter names to estimate OK
param_est_names <- c(param_proc_est_names, param_iv_est_names)
# all fixed parameters OK
param_fixed_names <- c(param_proc_fixed_names, param_vacc_fixed_names)
# all param names OK
param_names <- c(param_est_names, param_fixed_names)

# names of parameters that are rates (MAYBE) (because time 365 since timestep is year)
param_rates_in_days_names <- c("mu", "alpha", "gammaI", "gammaA", "rhoA") #muB

# Measurement model  -------------------------------------------------------

## density

## NegBinomial density (if k -> inf then becomes Poisson)
dmeas <- Csnippet("
  double mean_cases = epsilon * C;
  if (t > 2018)
    mean_cases = mean_cases * cas_def;
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
  if (t > 2018)
    mean_cases = mean_cases * cas_def;
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
cases_at_t_start <- cases %>% filter(dateToYears(date) <= t_start)
cases_at_t_start.string <- foreach(r = iter(cases_at_t_start, by = "row"),
                                   .combine = c) %do% {
                                     sprintf(" {%f, %f} ", r$time, r$cases)
                                   } %>% 
  str_c(collapse = ", \n")

matrix_cases_at_t_start.string <- str_c(sprintf("double cases_at_t_start[%i][%i] = {\n", nrow(cases_at_t_start), 2),
                                        cases_at_t_start.string,
                                        " \n };")


# Cases from other departement as a mobility rational
cases_other.string <- foreach(r = iter(cases_other_dept, by = "row"),
                                   .combine = c) %do% {
                                     sprintf(" {%f, %f} ", r$time, r$cases)
                                   } %>% 
  str_c(collapse = ", \n")

matrix_cases_other.string <- str_c(sprintf("double cases_other[%i][%i] = {\n", nrow(cases_other_dept), 2),
                                        cases_other.string,
                                        " \n };")



initalizeStates <- Csnippet("
  A     = nearbyint((1-sigma)/sigma  * 1/epsilon * cases_at_t_start[n_cases_start-1][1]/7 * 365 /(mu+gammaA));
  I     = nearbyint(1/epsilon * cases_at_t_start[n_cases_start-1][1]/7 * 365 /(mu+alpha+gammaI))  ;  // Steady state, DP says its correct.
  double R0[2] = {0,0};
  //FILE *fp;
  //fp = fopen('Ouuuutput.txt', 'w');

  double B_acc = 0;
  double rhoI = rhoA * XrhoI;
  double thetaA = thetaI * XthetaA;
  for(int i = 0; i < n_cases_start; i++){
    R0[0] +=                   cases_at_t_start[i][1]/epsilon  * exp((cases_at_t_start[i][0] - t_start)  * (rhoI+mu)); /* because t_i in past so t_ - t_0 negative */
    R0[1] += (1-sigma)/sigma * cases_at_t_start[i][1]/epsilon  * exp((cases_at_t_start[i][0] - t_start)  * (rhoA+mu));
    B_acc += (thetaA * (1-sigma)/sigma * cases_at_t_start[i][1]/epsilon + thetaI * cases_at_t_start[i][1]/epsilon) *
              (1 + lambdaR * pow(0.024, r)) * D * exp((cases_at_t_start[i][0] - t_start)  * mu_B);
    
    //fprintf(fp, '%f %f %f', cases_at_t_start[i][0] - t_start, cases_at_t_start[i][0], t_start);
  }
 //fclose(fp);

  B = B_acc;
  RI1   = nearbyint(R0[0]/3);
  RI2   = nearbyint(R0[0]/3);
  RI3   = nearbyint(R0[0]/3);
  RA1   = nearbyint(R0[1]/3);
  RA2   = nearbyint(R0[1]/3);
  RA3   = nearbyint(R0[1]/3);
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
  B   = (I * thetaI/mu_B + A * thetaA/mu_B) * D * (1 + lambdaR * pow(0.024, r)); // TODO custom initial conditions equivalent to the 'forcing' in the continous model
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
  Tcas_def = logit(cas_def);
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
  Tcas_def = expit(cas_def);
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
param_fixed[param_proc_fixed_names] <- as.numeric(param_proc_fixed)
param_fixed['cases_ext'] <- 1

# Initialize the parameters to estimate (just initial guesses)
param_est <- set_names(seq_along(param_est_names) * 0, param_est_names)
param_est["sigma"] <- .25               # Fixed
param_est["rhoA"] <- 1/(365*8)*365.25   # Fixed
param_est["XrhoI"] <- 1                 # Fixed
param_est["gammaA"] <- 182.625          # Fixed 	
param_est["gammaI"] <- 182.625          # Fixed
param_est["Rtot_0"] <- 0.35             # Useless

# Estimated parameters for all dept:
param_est["betaB"] <- 0.296517
param_est["foi_add"] <- 0.0
# Estimated parameters for Artibonite:
param_est["mu_B"] <-  159.593008
param_est["XthetaA"] <- 0.033616
param_est["thetaI"] <- 0.000767
param_est["lambdaR"] <- 0.902726
param_est["r"] <- 1.041372
param_est["std_W"] <- 0.007852
param_est["epsilon"] <- 0.975316
param_est["k"] <- 379.266151
param_est["cas_def"] <- 0.5



cases_ext_mean <- cases_other_dept %>% filter(time > t_start)
cases_ext_mean <- mean(cases_ext_mean$cases)

# rate of simulation in fractions of years
dt_yrs <- 1/365.25 * .2

# adjust the rate parameters depending on the integration delta time in years (some parameter inputs given in days)
params <- c(param_est, param_fixed)  
#params[param_rates_in_days_names] <- params[param_rates_in_days_names] * 365.25
# Was mu", "alpha", "gammaI", "gammaA", "rhoA"
rain <- rain  %>% 
    filter(time > (t_start - 0.01) & time < (t_end + 0.01)) %>% 
    select(time, rain_std) %>% 
    rename(rain = rain_std)

cases_covar <- cases_covar %>% 
    filter(time > (t_start - 0.01) & time < (t_end + 0.01)) %>% 
    select(time, cases) %>% 
    rename(cases_covar_c = cases)

covar <- full_join(rain, cases_covar)


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
  covar = covar,
  tcovar = "time",
  # names of state variables
  statenames = state_names,
  # names of accumulator variables to be re-initalized at each observation timestep 
  # (C for cases, W for the white noise just for plotting)
  zeronames = c("C", "W"),
  # names of paramters
  paramnames = param_names,
  # names of covariates
 # covarnames = "rain",
  initializer = initalizeStates,
  toEstimationScale = toEstimationScale,
  fromEstimationScale = fromEstimationScale,
  # global C definitions
  globals = str_c(
    sprintf("double t_vacc_start_alt = %f; double t_vacc_end_alt = %f;", t_vacc_start_alt, t_vacc_end_alt),
    sprintf("double r_v_alt = %f;", r_v_alt_year),
    sprintf("double p1d_alt = %f;", p1d_alt),
    sprintf("int n_cases_start = %i;",  nrow(cases_at_t_start)),
    sprintf("int n_cases_other = %i;",  nrow(cases_other_dept)),
    sprintf("double t_start = %f;",  t_start),
    sprintf("double t_end = %f;",  t_end),
    derivativeBacteria.c,
    matrix_cases_at_t_start.string,
    matrix_cases_other.string,
    eff_v.c,
    sep = " ")
)

# save pomp object for further use
save(sirb_cholera, file = paste0(output_dir, departement, "/sirb_cholera_pomped_", departement, ".rda"))
