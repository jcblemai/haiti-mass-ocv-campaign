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

# scenario     <- "S0"
# nsim         <- 100
# t_vacc_start <- list()
# t_vacc_end <- list()
# p1d_reg<- list()
# r_v_year <- list()
# t_vacc_start$Artibonite <- "2024-01-11"
# t_vacc_end$Artibonite   <- "2029-01-09"
# p1d_reg$Artibonite     <- 1.0
# r_v_year$Artibonite     <- 0.0003457414265972485
# cases_ext    <- 1
# t_vacc_start$Centre <- "2019-01-12"
# t_vacc_end$Centre   <- "2024-01-11"
# p1d_reg$Centre     <- 1.0
# r_v_year$Centre     <- 0.00014934941524298614
# cases_ext    <- 1
# t_vacc_start$Grande_Anse <- "2064-01-01"
# t_vacc_end$Grande_Anse   <- "2068-12-30"
# p1d_reg$Grande_Anse     <- 1.0
# r_v_year$Grande_Anse     <- 9.372434525767405e-05
# cases_ext    <- 1
# t_vacc_start$Nippes <- "2049-01-04"
# t_vacc_end$Nippes   <- "2054-01-03"
# p1d_reg$Nippes     <- 1.0
# r_v_year$Nippes     <- 6.855191716307418e-05
# cases_ext    <- 1
# t_vacc_start$Nord <- "2039-01-07"
# t_vacc_end$Nord   <- "2044-ß01-06"
# p1d_reg$Nord     <- 1.0
# r_v_year$Nord     <- 0.0002135815759501876
# cases_ext    <- 1
# t_vacc_start$Nord_Est <- "2054-01-03"
# t_vacc_end$Nord_Est   <- "2059-01-02"
# p1d_reg$Nord_Est     <- 1.0
# r_v_year$Nord_Est     <- 7.884736340116735e-05
# cases_ext    <- 1
# t_vacc_start$Nord_Ouest <- "2034-01-08"
# t_vacc_end$Nord_Ouest   <- "2039-01-07"
# p1d_reg$Nord_Ouest     <- 1.0
# r_v_year$Nord_Ouest     <- 0.0001458612279158269
# cases_ext    <- 1
# t_vacc_start$Ouest <- "2029-01-09"
# t_vacc_end$Ouest   <- "2034-01-08"
# p1d_reg$Ouest     <- 1.0
# r_v_year$Ouest     <- 0.0008064929665035424
# cases_ext    <- 1
# t_vacc_start$Sud <- "2044-01-06"
# t_vacc_end$Sud   <- "2049-01-04"
# p1d_reg$Sud     <- 1.0
# r_v_year$Sud     <- 0.00015510135188780548
# cases_ext    <- 1
# t_vacc_start$Sud_Est <- "2059-01-02"
# t_vacc_end$Sud_Est   <- "2064-01-01"
# p1d_reg$Sud_Est     <- 1.0
# r_v_year$Sud_Est     <- 0.00012660685015481463
# cases_ext    <- 1



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



# liks_stoch <- read_csv("output_17-17-mobinf/Haiti_OCVAll-param_logliks-10-l2.csv") 
# liks <- liks_stoch
# best_param <- liks %>% 
#   arrange(desc(loglik)) %>% 
#   slice(1)  %>% 
#   arrange(desc(loglik)) %>% 
#   ungroup %>% 
#   left_join(liks_stoch)
load(paste0(output_dir, "/sirb_cholera_pomped_all.rda"))
#params <- unlist(best_param[names(coef(sirb_cholera))]) %>% as.double()
#names(params) <- names(best_param[names(coef(sirb_cholera))])
# pomp::coef(sirb_cholera) <- params   # UNCOMMENT TO CHANGE TO CALIBRATION PARAMETERS

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


all_rain <- read_csv("haiti-data/fromAzman/rainfall.csv") %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date)) %>%
  filter(time > t_start - 0.01 & time < (t_end + 0.01))


for (dp in departements) {
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

}

all_rain_forecast <- read_csv("haiti-data/proj/rainfall.csv") %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date)) %>%
  filter(time > t_start - 0.01 & time < (t_forecast + 0.01))

for (dp in departements) {
  rain_forecast <- read_csv("haiti-data/proj/rainfall.csv")  %>% 
    gather(dep, rain, -date) %>% 
    group_by(dep) %>% 
    ungroup() %>% 
    filter(dep == dp) %>%
    mutate(date = as.Date(date, format = "%Y-%m-%d"),
           time = dateToYears(date)) %>%
    filter(time > t_start - 0.01 & time < (t_forecast + 0.01)) %>%
    mutate(max_rain = all_rain[paste0('max_rain',gsub('-','_', dp))][[1]][1],
           rain_std = rain/all_rain[paste0('max_rain',gsub('-','_', dp))][[1]][1])
  
  all_rain_forecast <- cbind(all_rain_forecast, placeholder2 = rain_forecast$rain_std)
  all_rain_forecast <- cbind(all_rain_forecast, placeholder = rain_forecast$max_rain)
  names(all_rain_forecast)[names(all_rain_forecast) == "placeholder2"] <- paste0('rain_std', gsub('-','_',dp))
  names(all_rain_forecast)[names(all_rain_forecast) == "placeholder"] <- paste0('max_rain', gsub('-','_',dp))

}

time_forecast <- dateToYears(seq.Date(yearsToDate(t_start), yearsToDate(t_forecast), by = "1 week"))

doMC::registerDoMC(8)

# function to simulate from a given set of paramters
simulatePOMP <- function(nsim, seed = 199919L) {
  for (dp in departements) {
    coef(sirb_cholera)[paste0("t_vacc_start", gsub('-','_', dp))] <- dateToYears(as.Date(t_vacc_start[gsub('-','_', dp)][[1]][1]))
    coef(sirb_cholera)[paste0("t_vacc_end", gsub('-','_', dp))] <- dateToYears(as.Date(t_vacc_end[gsub('-','_', dp)][[1]][1]))
    coef(sirb_cholera)[paste0("p1d_reg", gsub('-','_', dp))] <-  as.numeric(p1d_reg[gsub('-','_', dp)][[1]][1])
    coef(sirb_cholera)[paste0("r_v_year", gsub('-','_', dp))] <- as.numeric(r_v_year[gsub('-','_', dp)][[1]][1])
  }
  coef(sirb_cholera)["cases_ext"] <- as.numeric(cases_ext)

  pomp::simulate(sirb_cholera, nsim = nsim, as.data.frame = T , include.data = TRUE, seed = runif(1,1,10000), times = time_forecast) -> projec
  save(projec, file = sprintf("output/Simulations/Haiti_OCV_Projection-allDep-%i-%s.rda", nsim, scenario))
  
  projec %>%
    as_tibble() %>% 
    mutate(isdata = sim == "data") %>%
    gather(variable, value, -time, -sim, -isdata) %>% 
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
                     covar = all_rain_forecast,
                     tcovar = "time")

# run simulations for each model
sim_stochastic <- simulatePOMP(nsim = nsim)


