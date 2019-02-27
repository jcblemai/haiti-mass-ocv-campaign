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

# TODO comment to run in python
#output_dir <- "output_14-21-newdata/"
#departement <- 'Artibonite'
#run_level <- 3
#nsim <- 1000
#cases_ext <- 1 
#t_vacc_start <- '2010-01-01'
#t_vacc_end  <- '2010-01-01'
#p1d_reg <- 0 
#r_v_year <- 0

liks_stoch <- read_csv(sprintf("%s%s/Haiti_OCV-%s-param_logliks-10-l%i.csv",output_dir, departement, departement, run_level)) 
liks <- liks_stoch

# get MLE paramter sets 
best_param <- liks %>% 
  arrange(desc(loglik)) %>% 
  slice(1)  %>% 
  arrange(desc(loglik)) %>% 
  ungroup %>% 
  left_join(liks_stoch)

load(paste0(output_dir, departement, "/sirb_cholera_pomped_", departement, ".rda"))

params <- unlist(best_param[names(coef(sirb_cholera))]) %>% as.double()
names(params) <- names(best_param[names(coef(sirb_cholera))])


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


#case_dates <- cases %>%  filter(time > t_start - 0.01)

# check if weekly reports are missing
#missing_dates <- setdiff(seq.Date(yearsToDate(t_start), yearsToDate(t_forecast), by = "1 week"), case_dates$date) %>% as.Date(origin = as.Date("1970-01-01"))

#if (length(missing_dates) > 0) {
# fill in the data
#  cases %<>% 
#    bind_rows(slice(cases, 1:2) %>% 
#                mutate(date = missing_dates,
#                       cases = NA,
#                       time = dateToYears(date))) %>%
#    arrange(time)
#}

# Compare outputs ---------------------------------------------------------

doMC::registerDoMC(8)

# function to simulate from a given set of paramters
simulatePOMP <- function(params, nsim, seed = 199919L) {
  pomp::coef(sirb_cholera) <- params
  coef(sirb_cholera)["t_vacc_start"] <- dateToYears(as.Date(t_vacc_start))
  coef(sirb_cholera)["t_vacc_end"] <- dateToYears(as.Date(t_vacc_end))
  coef(sirb_cholera)["p1d_reg"] <-  as.numeric(p1d_reg)
  coef(sirb_cholera)["r_v_year"] <- as.numeric(r_v_year)
  coef(sirb_cholera)["cases_ext"] <- as.numeric(cases_ext)

  pomp::simulate(sirb_cholera, nsim = nsim, as.data.frame = T , include.data = TRUE, seed = runif(1,1,10000), times = time_forecast) -> calib
  save(calib, file = "calib.rda")
  
  calib %>%
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

# Build new covariate:
cases_covar <- read_csv("covar_mob.csv")  %>% 
  gather(dep, cases, -date) %>% 
  filter(dep != departement) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date))

cases_covar <- aggregate(cases_covar$cases, by=list(Category=cases_covar$time), FUN=sum, na.rm=TRUE, na.action=NULL) %>%
  mutate(time = Category) %>%
  mutate(cases = x)

rain <- rain_forecast  %>% 
  filter(time > (t_start - 0.01) & time < (t_forecast + 0.01)) %>% 
  select(time, rain_std) %>% 
  rename(rain = rain_std)

cases_covar <- cases_covar %>% 
  filter(time > (t_start - 0.01) & time < (t_forecast + 0.01)) %>% 
  select(time, cases) %>% 
  rename(cases_covar_c = cases)

covar <- full_join(rain, cases_covar) %>% # %>% #default direction down
  fill(cases_covar_c, .direction = "up") %>% fill(cases_covar_c)

sirb_cholera <- pomp(sirb_cholera,
                     covar = covar,
                     tcovar = "time")

# run simulations for each model
sim_stochastic <- simulatePOMP(params, nsim = nsim)

# tidy tibble for merger
sim_stochastic_quantiles <- sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "cases", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)# %>% 
#bind_cols(sim_stochastic %>%                                         #Add a column with the data
#            filter(variable == "cases", isdata == "data") %>% 
#            select(mean) %>% 
#              rename(cases = mean))


# create objects for python
cases <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "cases", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

A <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "A", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

B <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "B", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

C <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "C", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

I <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "I", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

S <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "S", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

W <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "W", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)


RA1 <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "RA1", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)
RA2 <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "RA2", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

RA3 <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "RA3", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

RI1 <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "RI1", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

RI2 <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "RI2", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

RI3 <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "RI3", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)


VSd <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VSd", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRI1d <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRI1d", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRI2d <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRI2d", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRI3d <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRI3d", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRA1d <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRA1d", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRA2d <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRA2d", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRA3d <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRA3d", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)


# double doses:

VSdd <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VSdd", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRI1dd <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRI1dd", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRI2dd <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRI2dd", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRI3dd <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRI3dd", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRA1dd <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRA1dd", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRA2dd <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRA2dd", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRA3dd <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRA3dd", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

# PLOT BOTH
# simcol <- "#175CD6"
# datacol <- "#ED0000"
# 
# p.sim <- ggplot(data = sim_stochastic_quantiles,
#                 aes(x = date))+
#   geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.1, color = simcol, fill = simcol) +
#   geom_line(aes(y = q50), color = simcol) +
#   geom_line(aes(y = mean), linetype = 2, color = simcol) +
#   #geom_line(aes(y = cases), color = datacol, lwd = 0.2) +
#   #geom_point(aes(y = cases), color = datacol, size = 0.8) +
#   #geom_text(data = psim_labels, aes (y = y, label = label), size = 7) +
#   #facet_grid(model~type) +
#   scale_x_date(date_labels = "%b-%y", expand = c(0,0), limits = as.Date(c(yearsToDate(t_start), yearsToDate(t_forecast)))) +
#   scale_y_continuous(expand = c(0,0))+
#   labs(y = "daily cholera cases", x = "date") +
#   theme(panel.grid.major = element_line(color = "lightgray"),
#         panel.background = element_rect(fill = "white"),
#         axis.line = element_line(color = "black"),
#         strip.text = element_blank(),
#         axis.title = element_text())
# 
# print(p.sim)
# 
# sim_stochastic_quantiles_all <- sim_stochastic %>% 
#   mutate(date = as.Date(round_date(date))) %>% 
#   filter(variable == "cases" |  isdata == "simulation") %>% 
#   select(-isdata) 
# 
# p.all <- ggplot(data = sim_stochastic_quantiles_all,
#                 aes(x = date))+
#   geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.1, color = simcol, fill = simcol) +
#   geom_line(aes(y = q50), color = simcol) +
#   geom_line(aes(y = mean), linetype = 2, color = simcol) +
#   facet_wrap(~variable, scales = "free_y") 
# scale_x_date(date_labels = "%b-%y", expand = c(0,0), limits = as.Date(c("2014-03-01", "2018-07-14")))

#print(p.all)




VSd_alt <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VSd_alt", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRI1d_alt <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRI1d_alt", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRI2d_alt <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRI2d_alt", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRI3d_alt <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRI3d_alt", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRA1d_alt <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRA1d_alt", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRA2d_alt <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRA2d_alt", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRA3d_alt <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRA3d_alt", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)


# double doses:

VSdd_alt <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VSdd_alt", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRI1dd_alt <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRI1dd_alt", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRI2dd_alt <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRI2dd_alt", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRI3dd_alt <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRI3dd_alt", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRA1dd_alt <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRA1dd_alt", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRA2dd_alt <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRA2dd_alt", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)

VRA3dd_alt <-sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "VRA3dd_alt", isdata == "simulation") %>% 
  select(-isdata, -time, -variable)


