library(tidyverse)
library(magrittr)
library(ggthemes)
library(GGally)
library(foreach)
library(itertools)
library(pomp)
library(lubridate)
Sys.setlocale("LC_ALL","C")
Sys.setenv(TZ='GMT')

# STOCHASTIC MODEL
# Auxillary functions
dateToYears <- function(date, origin = as.Date("2014-01-01"), yr_offset = 2014) {
  julian(date, origin = origin)/365.25 + yr_offset
}

yearsToDate <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
  as.Date((year_frac - yr_offset) * 365.25, origin = origin)
}

yearsToDateTime <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
  as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
}
output_dir <- "output/"

load(paste0(output_dir, "/sirb_cholera_pomped_all.rda"))

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




# input parameters to the model
input_parameters <- yaml::read_yaml("haiti-data/input_parameters.yaml")
# Start and end dates of epidemic
t_start <- dateToYears(as.Date(input_parameters$t_start))
t_end <- dateToYears(as.Date(input_parameters$t_end))

t_forecast <- dateToYears(as.Date("2029-12-21"))


rain4max <- read_csv("haiti-data/fromAzman/rainfall.csv")  %>% 
  gather(dep, rain, -date) %>% 
  group_by(dep) %>% 
  ungroup() %>% 
  filter(dep == departement) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date)) %>%
  filter(time > t_start - 0.01 & time < (t_end + 0.01)) %>%
  mutate(max_rain = max(rain), rain_std = rain/max_rain) 


max_rain<- rain4max %>% select(max_rain) %>% max()


rain_forecast <- read_csv("haiti-data/proj/rainfall.csv")  %>% 
  gather(dep, rain, -date) %>% 
  group_by(dep) %>% 
  ungroup() %>% 
  filter(dep == departement) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date)) %>%
  filter(time > t_start - 0.01 & time < (t_forecast + 0.01)) %>%
  mutate(max_rain = max_rain, rain_std = rain/max_rain) 


time_forecast <- dateToYears(seq.Date(yearsToDate(t_start), yearsToDate(t_forecast), by = "1 week"))

doMC::registerDoMC(8)

# function to simulate from a given set of paramters
simulatePOMP <- function(nsim, seed = 199919L) {
  for (dp in departements) {
    coef(sirb_cholera)[paste0("t_vacc_start", gsub('-','_', dp))] <- dateToYears(as.Date(t_vacc_start$dp))
    coef(sirb_cholera)[paste0("t_vacc_end", gsub('-','_', dp))] <- dateToYears(as.Date(t_vacc_end$dp))
    coef(sirb_cholera)[paste0("p1d_reg", gsub('-','_', dp))] <-  as.numeric(p1d_reg$dp)
    coef(sirb_cholera)[paste0("r_v_year", gsub('-','_', dp))] <- as.numeric(r_v_year$dp)
  }
  coef(sirb_cholera)["cases_ext"] <- as.numeric(cases_ext)
  coef(sirb_cholera)["cas_def"] <- 0.36
  
  pomp::simulate(sirb_cholera, nsim = nsim, as.data.frame = T , include.data = TRUE, seed = runif(1,1,10000), times = time_forecast) -> projec
  save(projec, file = sprintf("output/Simulations/Haiti_OCV_Projection-allDep-%i-%s.rda", nsim, scenario))
  
  projec %>%
    as_tibble() %>% 
    mutate(isdata = sim == "data") %>%
    gather(variable, value, -time, -rain, -sim, -isdata) %>% 
    group_by(time, isdata, variable) %>% 
    summarise( q05 = quantile(value, 0.025, na.rm = T),
               mean = mean(value, na.rm = T),
               q50 = quantile(value, 0.5, na.rm = T),
               q95 = quantile(value, 0.975, na.rm = T)) %>% 
    ungroup %>% 
    mutate(isdata = ifelse(isdata, "data", "simulation"),
           date = yearsToDateTime(time)) %>% 
    filter(date >= yearsToDate(sirb_cholera@t0))
  
}

sirb_cholera <- pomp(sirb_cholera,
                     covar = rain,
                     tcovar = "time")

# run simulations for each model
sim_stochastic <- simulatePOMP(nsim = nsim)

# tidy tibble for merger
sim_stochastic_quantiles <- sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "cases", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)# %>% 

