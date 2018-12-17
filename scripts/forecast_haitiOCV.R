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

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  # default departement
  args[1] = "Artibonite"
  args[2] = 3
} else if (length(args)==1) {
  args[2] = 1
}

departement <- args[1]
run_level <- as.integer(args[2])
nsim = 10


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

params <-unlist(best_param[names(coef(sirb_cholera))])









# input parameters to the model
input_parameters <- yaml::read_yaml("haiti-data/input_parameters.yaml")
# Start and end dates of epidemic
t_start <- dateToYears(as.Date(input_parameters$t_start))
t_end <- dateToYears(as.Date(input_parameters$t_end))

t_sf <- dateToYears(as.Date(input_parameters$t_end))


t_forecast <- dateToYears(as.Date("2029-12-20"))


rain_forecast <- read_csv("haiti-data/proj/rainfall.csv")  %>% 
  gather(dep, rain, -date) %>% 
  group_by(dep) %>% 
  ungroup() %>% 
  filter(dep == departement) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date)) %>%
  filter(time > t_start - 0.01 & time < (t_forecast + 0.01)) %>%
  mutate(max_rain = max(rain), rain_std = rain/max_rain) 


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
  #pomp::coef(sirb_cholera_forecast) <- params
  
  pomp::simulate(sirb_cholera, nsim = nsim, as.data.frame = T , include.data = TRUE, seed = seed, times = time_forecast) -> calib
  #pomp::simulate(sirb_cholera_forecast, nsim = nsim, as.data.frame = T , include.data = TRUE, seed = seed, times = time_forecast) -> proj
  
  #rbind(calib, proj)%>% 
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



# sirb_cholera_maxrain <- pomp(sirb_cholera,
#                              covar = rain %>%
#                                filter(time >= min(sirb_cholera@tcovar) & time <= max(sirb_cholera@tcovar)) %>% 
#                                select(time, rain_std) %>%
#                                rename(rain = rain_std),
#                              tcovar = "time")


# New ofject with new horizon:
#sirb_cholera_forecast <- sirb_cholera
#time(sirb_cholera_forecast) <- time_forecast
#timezero(sirb_cholera_forecast) <- max(time(sirb_cholera))
#covar(sirb_cholera_forecast) <- rain_forecast
#sirb_cholera_forecast <- pomp(sirb_cholera_forecast,
#                              covar = rain_forecast,
#                              tcovar = "time")

sirb_cholera <- pomp(sirb_cholera,
                              covar = rain_forecast,
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


# create ojects for python

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
# #print(p.sim)
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


