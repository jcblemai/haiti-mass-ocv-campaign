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
Sys.setlocale("LC_ALL", "C")
output_dir <- "output/"

# function to convert dates to fractions of years for model
dateToYears <-
  function(date,
           origin = as.Date("2014-01-01"),
           yr_offset = 2014) {
    julian(date, origin = origin) / 365.25 + yr_offset
  }

yearsToDate <-
  function(year_frac,
           origin = as.Date("2014-01-01"),
           yr_offset = 2014.0) {
    as.Date((year_frac - yr_offset) * 365.25, origin = origin)
  }

yearsToDateTime <-
  function(year_frac,
           origin = as.Date("2014-01-01"),
           yr_offset = 2014.0) {
    as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
  }

departements <-
  c(
    'Artibonite',
    'Centre',
    'Grande_Anse',
    'Nippes',
    'Nord',
    'Nord-Est',
    'Nord-Ouest',
    'Ouest',
    'Sud',
    'Sud-Est'
  )

# define stat variable names OK
state_names <-
  c(
    "S",
    "I",
    "A",
    "RI1",
    "RI2",
    "RI3",
    "RA1",
    "RA2",
    "RA3",
    "VSd",
    "VRI1d",
    "VRI2d",
    "VRI3d",
    "VRA1d",
    "VRA2d",
    "VRA3d",
    "VSdd",
    "VRI1dd",
    "VRI2dd",
    "VRI3dd",
    "VRA1dd",
    "VRA2dd",
    "VRA3dd",
    "VSd_alt",
    "VRI1d_alt",
    "VRI2d_alt",
    "VRI3d_alt",
    "VRA1d_alt",
    "VRA2d_alt",
    "VRA3d_alt",
    "VSdd_alt",
    "VRI1dd_alt",
    "VRI2dd_alt",
    "VRI3dd_alt",
    "VRA1dd_alt",
    "VRA2dd_alt",
    "VRA3dd_alt",
    "B",
    "C",
    "W"
  )


params_common <-
  c(
    "sigma",
    "mu_B",
    "thetaI",
    "XthetaA",
    "lambdaR",
    "r",
    "gammaI",
    "gammaA",
    "rhoA",
    "XrhoI",
    "epsilon",
    "k",
    "std_W",
    "cas_def",
   "Rtot_0",
    "mu",
   "alpha", 
   "cases_ext")

## fixed process model parameters  OK
params_diff <- c("foi_add", "betaB", "H", "D","t_vacc_start", "t_vacc_end", "p1d_reg", "r_v_year",
    "t_vacc_start_alt", "t_vacc_end_alt", "p1d_reg_alt", "r_v_year_alt")

# input parameters to the model
input_parameters <-
  yaml::read_yaml("haiti-data/input_parameters.yaml")

# Start and end dates of epidemic
t_start <- dateToYears(as.Date(input_parameters$t_start))
t_end   <- dateToYears(as.Date(input_parameters$t_end))

all_state_names <- c('IncidenceAll', 'DosesAll', 'CasesAll')
all_param_names <- params_common
all_matrix_cases_at_t_start.string <- ""
all_matrix_cases_other.string <- ""
all_rain <- read_csv("haiti-data/fromAzman/rainfall.csv") %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date)) %>%
  filter(time > t_start - 0.01 & time < (t_end + 0.01))

all_cases <- read_csv("haiti-data/fromAzman/cases_corrected.csv")  %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date))

for (dp in departements) {
  cases <- read_csv("haiti-data/fromAzman/cases_corrected.csv")  %>%
    gather(dep, cases,-date) %>%
    filter(dep == dp) %>%
    mutate(date = as.Date(date, format = "%Y-%m-%d"),
           time = dateToYears(date))

  
  rain <- read_csv("haiti-data/fromAzman/rainfall.csv")  %>%
    gather(dep, rain,-date) %>%
    group_by(dep) %>%
    ungroup() %>%
    filter(dep == dp) %>%
    mutate(date = as.Date(date, format = "%Y-%m-%d"),
           time = dateToYears(date)) %>%
    filter(time > t_start - 0.01 & time < (t_end + 0.01)) %>%
    mutate(max_rain = max(rain), rain_std = rain / max_rain)

  all_rain <- cbind(all_rain, placeholder = rain$max_rain)
  all_rain <- cbind(all_rain, placeholder2 = rain$rain_std)
  names(all_rain)[names(all_rain) == "placeholder"] <- paste0('max_rain', gsub('-','_',dp))
  names(all_rain)[names(all_rain) == "placeholder2"] <- paste0('rain_std', gsub('-','_',dp))
  
  all_cases <- cbind(all_cases, placeholder = cases$cases)
  names(all_cases)[names(all_cases) == "placeholder"] <- paste0('cases', gsub('-','_',dp))
  
  
  cases_other_dept <-
    read_csv("haiti-data/fromAzman/cases_corrected.csv")  %>%
    gather(dep, cases,-date) %>%
    filter(dep != dp) %>%
    mutate(date = as.Date(date, format = "%Y-%m-%d"),
           time = dateToYears(date))
  
  cases_other_dept <-
    aggregate(
      cases_other_dept$cases,
      by = list(Category = cases_other_dept$time),
      FUN = sum,
      na.rm = TRUE,
      na.action = NULL
    ) %>%
    mutate(time = Category) %>%
    mutate(cases = x)
  
  cases_at_t_start <- cases %>% filter(dateToYears(date) <= t_start)
  cases_at_t_start.string <-
    foreach(r = iter(cases_at_t_start, by = "row"),
            .combine = c) %do% {
              sprintf(" {%f, %f} ", r$time, r$cases)
            } %>%
    str_c(collapse = ", \n")
  
  matrix_cases_at_t_start.string <-
    str_c(
      sprintf(
        "double cases_at_t_start%s[%i][%i] = {\n",
        gsub('-', '_', dp),
        nrow(cases_at_t_start),
        2
      ),
      cases_at_t_start.string,
      " \n };"
    )
  
  # Cases from other departement as a mobility rational
  cases_other.string <-
    foreach(r = iter(cases_other_dept, by = "row"),
            .combine = c) %do% {
              sprintf(" {%f, %f} ", r$time, r$cases)
            } %>%
    str_c(collapse = ", \n")
  
  matrix_cases_other.string <-
    str_c(sprintf("double cases_other%s[%i][%i] = {\n", gsub('-', '_', dp), nrow(cases_other_dept), 2),
          cases_other.string,
          " \n };")
  
  all_matrix_cases_at_t_start.string = str_c(all_matrix_cases_at_t_start.string, matrix_cases_at_t_start.string)
  all_matrix_cases_other.string = str_c(all_matrix_cases_other.string, matrix_cases_other.string)


  all_state_names = append(all_state_names, lapply(state_names, paste0, gsub('-', '_',dp)))
  all_param_names = append(all_param_names, lapply(params_diff, paste0, gsub('-', '_',dp)))
  
}
all_state_names = unlist(all_state_names)
all_param_names = unlist(all_param_names)

all_params <-
  set_names(seq_along(all_param_names) * 0, all_param_names)

all_params["sigma"] <- .25               # Fixed
all_params["rhoA"] <- 1 / (365 * 8) * 365.25   # Fixed
all_params["XrhoI"] <- 1                 # Fixed
all_params["gammaA"] <- 182.625          # Fixed
all_params["gammaI"] <- 182.625          # Fixed
all_params["Rtot_0"] <- 0.35             # Useless


all_params['cases_ext'] <- 1
all_params['mu'] <-  0.01586625546  
all_params['alpha'] <- 1.461

# From calibration
all_params["mu_B"] <-  133.19716102404308
all_params["XthetaA"] <- 0.0436160721505241
all_params["thetaI"] <- 3.4476623459780395e-4
all_params["lambdaR"] <- 0.2774237712085347
all_params["r"] <- 0.31360358752214235
all_params["std_W"] <- 0.008172280355938182
all_params["epsilon"] <- 0.9750270707877388
all_params["k"] <- 101.2215999283583
all_params["cas_def"] <- 0.10

all_params["betaBArtibonite"] =   0.516191 
all_params["betaBSud_Est"] =      1.384372 
all_params["betaBNippes"] =       2.999928 
all_params["betaBNord_Est"] =     3.248645 
all_params["betaBOuest"] =        0.090937 
all_params["betaBCentre"] =       1.977686 
all_params["betaBNord"] =         0.589541 
all_params["betaBSud"] =          1.305966 
all_params["betaBNord_Ouest"] =   1.141691 
all_params["betaBGrande_Anse"] =  2.823539 

all_params["foi_addArtibonite"] =    1.530994e-06  
all_params["foi_addSud_Est"] =     6.105491e-07
all_params["foi_addNippes"] =       3.056857e-07
all_params["foi_addNord_Est"] =     8.209611e-07
all_params["foi_addOuest"] =       1.070717e-06 
all_params["foi_addCentre"] =      0.0000106504579266415
all_params["foi_addNord"] =         5.319736e-07
all_params["foi_addSud"] =          1.030357e-06   
all_params["foi_addNord_Ouest"] =   5.855759e-07
all_params["foi_addGrande_Anse"] =  8.762740e-07 


for (dp in departements) {
  populations  <- unlist(flatten(input_parameters["population"]))
  densities <- unlist(flatten(input_parameters["density"]))
  all_params[paste0('H', gsub('-','_',dp))] <- populations[dp]
  all_params[paste0('D', gsub('-','_',dp))] <- densities[dp]
  p1d_alt_year  <- unlist(flatten(input_parameters["p1d_alt_year"]))
  nb_doses_alt_year <- unlist(flatten(input_parameters["nb_doses_alt_year"]))
  t_vacc_start_alt  <- unlist(flatten(input_parameters["t_vacc_start_alt"]))
  t_vacc_end_alt <- unlist(flatten(input_parameters["t_vacc_end_alt"]))
  # Vaccination information:
  t_vacc_start_alt = dateToYears(as.Date(t_vacc_start_alt[dp]))
  t_vacc_end_alt   = dateToYears(as.Date(t_vacc_end_alt[dp]))
  r_v_alt_year = nb_doses_alt_year[dp] / (t_vacc_end_alt - t_vacc_start_alt)
  p1d_alt = p1d_alt_year[dp]
  
  all_params[paste0("t_vacc_start_alt", gsub('-','_',dp))] = t_vacc_start_alt
  all_params[paste0("t_vacc_end_alt",gsub('-','_',dp))] = t_vacc_end_alt
  all_params[paste0("p1d_reg_alt", gsub('-','_',dp))] = p1d_alt 
  all_params[paste0("r_v_year_alt", gsub('-','_',dp))] = r_v_alt_year

}

initalizeStatesTemplate =   "
A%s     = nearbyint((1-sigma)/sigma  * 1/epsilon * cases_at_t_start%s[n_cases_start-1][1]/7 * 365 /(mu+gammaA));
I%s     = nearbyint(1/epsilon * cases_at_t_start%s[n_cases_start-1][1]/7 * 365 /(mu+alpha+gammaI))  ;  // Steady state, DP says its correct.

R0[0] = 0;
R0[1] = 0;
B_acc = 0;

for(int i = 0; i < n_cases_start; i++){
R0[0] +=                   cases_at_t_start%s[i][1]/epsilon  * exp((cases_at_t_start%s[i][0] - t_start)  * (rhoI+mu)); /* because t_i in past so t_ - t_0 negative */
R0[1] += (1-sigma)/sigma * cases_at_t_start%s[i][1]/epsilon  * exp((cases_at_t_start%s[i][0] - t_start)  * (rhoA+mu));
B_acc += (thetaA * (1-sigma)/sigma * cases_at_t_start%s[i][1]/epsilon + thetaI * cases_at_t_start%s[i][1]/epsilon) *
(1 + lambdaR * pow(0.024, r)) * D%s * exp((cases_at_t_start%s[i][0] - t_start)  * mu_B);

}

B%s = B_acc;
RI1%s   = nearbyint(R0[0]/3);
RI2%s   = nearbyint(R0[0]/3);
RI3%s   = nearbyint(R0[0]/3);
RA1%s   = nearbyint(R0[1]/3);
RA2%s   = nearbyint(R0[1]/3);
RA3%s   = nearbyint(R0[1]/3);
if (A%s + I%s + RI1%s + RI2%s + RI3%s + RA1%s + RA2%s + RA3%s >= H%s)
{
  double R_tot = H%s - A%s - I%s - 100.0;
  if (R_tot <= 0)
  {
  I%s     = nearbyint(H%s - 100);
  A%s     = nearbyint(0);
  R_tot = nearbyint(0);
  }
  RI1%s   = nearbyint(sigma * R_tot/3.0);
  RI2%s   = nearbyint(sigma * R_tot/3.0);
  RI3%s   = nearbyint(sigma * R_tot/3.0);
  RA1%s   = nearbyint((1-sigma) * R_tot/3.0);
  RA2%s   = nearbyint((1-sigma) * R_tot/3.0);
  RA3%s   = nearbyint((1-sigma) * R_tot/3.0);
}
S%s   = nearbyint(H%s - A%s - I%s - RI1%s - RI2%s - RI3%s - RA1%s - RA2%s - RA3%s);
B%s   = (I%s * thetaI/mu_B + A%s * thetaA/mu_B) * D%s * (1 + lambdaR * pow(0.024, r)); // TODO custom initial conditions equivalent to the 'forcing' in the continous model
C%s   = 0;
W%s   = 0;
VSd%s = 0;
VRI1d%s = 0;
VRI2d%s = 0;
VRI3d%s = 0;
VRA1d%s = 0;
VRA2d%s = 0;
VRA3d%s = 0;
VSdd%s = 0;
VRI1dd%s = 0;
VRI2dd%s = 0;
VRI3dd%s = 0;
VRA1dd%s = 0;
VRA2dd%s = 0;
VRA3dd%s = 0;
VSd_alt%s = 0;
VRI1d_alt%s = 0;
VRI2d_alt%s = 0;
VRI3d_alt%s = 0;
VRA1d_alt%s = 0;
VRA2d_alt%s = 0;
VRA3d_alt%s = 0;
VSdd_alt%s = 0;
VRI1dd_alt%s = 0;
VRI2dd_alt%s = 0;
VRI3dd_alt%s = 0;
VRA1dd_alt%s = 0;
VRA2dd_alt%s = 0;
VRA3dd_alt%s= 0;


"

initalizeStatesAll = "double R0[2] = {0,0};
IncidenceAll = 0;
DosesAll = 0;
CasesAll = 0;
double B_acc = 0;
double rhoI = rhoA * XrhoI;
double thetaA = thetaI * XthetaA;"
for (dp in departements) {
  initalizeStatesAll = paste0(initalizeStatesAll, gsub('%s', gsub('-', '_', dp), initalizeStatesTemplate))
}

initalizeStates <- Csnippet(initalizeStatesAll)


# Build pomp object -------------------------------------------------------??
rmeasTemplate <-   "
  double mean_cases%s = epsilon * C%s;
  if (t > 2018)
  mean_cases%s = mean_cases%s * cas_def;
  cases%s = rnbinom_mu(k, mean_cases%s);
  "

rmeasAll = ""
for (dp in departements) {
  rmeasAll = paste0(rmeasAll, gsub('%s', gsub('-', '_', dp), rmeasTemplate))
}

rmeas <- Csnippet(rmeasAll)


## NegBinomial density (if k -> inf then becomes Poisson)
dmeasTemplate <- "
    double mean_cases%s = epsilon * C%s;
    if (t > 2018)
       mean_cases%s = mean_cases%s * cas_def;
    if (ISNA(cases%s)) {
       lik += (give_log) ? 0 : 1;
    } else {
       if (S%s < 10000) {
          lik += (give_log) ? -99999 : 1.0e-18;
       } else {
           lik += dnbinom_mu(cases%s, k, mean_cases%s, give_log) ;
       }
    }
"

dmeasAll = "lik = 0;"
for (dp in departements) {
  dmeasAll = paste0(dmeasAll, gsub('%s', gsub('-', '_', dp), dmeasTemplate))
}

dmeas <- Csnippet(dmeasAll)

# Process model ----------------------------------------------------------------- OK

sirb_file <- 'scripts/sirb_model_all_dept.c'
sirb_file_init <- 'scripts/sirb_model_all_dept_init.c'



sirb.rprocTemplate = readChar(sirb_file, file.info(sirb_file)$size)

sirb.rproc = readChar(sirb_file_init, file.info(sirb_file_init)$size)
for (dp in departements) {
  sirb.rproc = paste0(sirb.rproc, gsub('%s', gsub('-', '_', dp), sirb.rprocTemplate))
}

sirb.rproc = Csnippet(sirb.rproc)

# C function to compute the time-derivative of bacterial concentration OK
derivativeBacteria.c <- " double fB(int I, int A, double B,
double mu_B, double thetaI, double XthetaA, double lambdaR, double rain, double r, double D) {
double thetaA = thetaI * XthetaA;
double dB;
dB = -mu_B * B +  (1 + lambdaR * pow(rain, r)) * D * (thetaI * (double) I + thetaA * (double) A);
return(dB);
};
"
eff_v.c <- paste0(readChar('scripts/v_eff.c', file.info(sirb_file)$size), " double eff_v_1d(double t_since_vacc, int scenario) {
  if (t_since_vacc < 1) 
                  return eff_v_2d(t_since_vacc, scenario);
                  else
                  return 0;
                  };
                  ")

zeronameTemplate = c("C", "W")
zeronameAll = c('IncidenceAll', 'CasesAll')
for (dp in departements){
  zeronameAll = append(zeronameAll, lapply(zeronameTemplate, paste0, gsub('-', '_', dp)))
  
}
zeronameAll <- unlist(zeronameAll)


toEstimationScale <- "       Tmu_B = log(mu_B);
                              TthetaI = log(thetaI);
                              TXthetaA = logit(XthetaA);
                              TlambdaR = log(lambdaR);
                              Tr = log(r);
                              Tstd_W = log(std_W);
                              Tepsilon = logit(epsilon);
                              Tcas_def = logit(cas_def);
                              Tk = log(k);
                              "
toEstimationScaleTemplate <- "TbetaB%s = log(betaB%s);
Tfoi_add%s = log(foi_add%s);"

for (dp in departements) {
  toEstimationScale = paste0(toEstimationScale, gsub('%s', gsub('-', '_', dp), toEstimationScaleTemplate))
}
toEstimationScale <- Csnippet(toEstimationScale)

fromEstimationScale <- "      Tmu_B = exp(mu_B);
                                TthetaI = exp(thetaI);
                                TXthetaA = expit(XthetaA);
                                TlambdaR = exp(lambdaR);
                                Tr = exp(r);
                                Tstd_W = exp(std_W);
                                Tepsilon = expit(epsilon);
                                Tcas_def = expit(cas_def);
                                Tk = exp(k);
                                "
fromEstimationScaleTemplate <- "TbetaB%s = exp(betaB%s);
Tfoi_add%s = exp(foi_add%s);"

for (dp in departements) {
  fromEstimationScale = paste0(fromEstimationScale, gsub('%s', gsub('-', '_', dp), fromEstimationScaleTemplate))
}

fromEstimationScale <- Csnippet(fromEstimationScale)


# rate of simulation in fractions of years
dt_yrs <- 1 / 365.25 * .2
sirb_cholera <- pomp(
  # set data
  data = all_cases %>%
    filter(time > t_start & time < (t_end + 0.01)) %>% select(time,
                                                              casesArtibonite,
                                                              casesCentre,
                                                              casesGrande_Anse,
                                                              casesNippes,
                                                              casesNord,
                                                              casesNord_Est,
                                                              casesOuest,
                                                              casesSud,
                                                              casesSud_Est,
                                                              casesNord_Ouest)
  ,
  # time column
  times = "time",
  # initialization time
  t0 = t_start - dt_yrs,
  # paramter vector
  params = all_params,
  
  # process simulator
  rprocess = euler.sim(step.fun = sirb.rproc, delta.t = dt_yrs),
  # measurement model simulator
  rmeasure =  rmeas,
  # measurement model density
  dmeasure = dmeas,
  # covariates
  covar = all_rain,
  tcovar = "time",
  # names of state variables
  statenames = all_state_names,
  # names of accumulator variables to be re-initalized at each observation timestep
  # (C for cases, W for the white noise just for plotting)
  zeronames = zeronameAll,
  # names of paramters
  paramnames = all_param_names,
  # names of covariates
  initializer = initalizeStates,
  toEstimationScale = toEstimationScale,
  fromEstimationScale = fromEstimationScale,
  # global C definitions
  globals = str_c(
    sprintf("int n_cases_start = %i;",  nrow(cases_at_t_start)),
    sprintf("int n_cases_other = %i;",  nrow(cases_other_dept)),
    sprintf("double t_start = %f;",  t_start),
    sprintf("double t_end = %f;",  t_end),
    derivativeBacteria.c,
    all_matrix_cases_at_t_start.string,
    all_matrix_cases_other.string,
    eff_v.c,
    sep = " "
  )
)

# save pomp object for further use
save(
  sirb_cholera,
  file = paste0(
    output_dir,
    "/sirb_cholera_pomped_all.rda"
  )
)
