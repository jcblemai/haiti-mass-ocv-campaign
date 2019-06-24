library(tidyverse)
library(magrittr)
library(ggthemes)
library(GGally)
library(foreach)
library(itertools)
library(pomp)
library(zoo)

library(lubridate)
Sys.setlocale("LC_ALL", "C")
Sys.setenv(TZ = 'GMT')

# STOCHASTIC MODEL
# Auxillary functions
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

team <- 'EPFL'
folder = '2019-05-27 Delivrable/Simulations/'
#folder = 'output/Simulations/'
load("output/sirb_cholera_pomped_all.rda")


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

df <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(df) <- departements

#all_data <- data.frame(departements)
nsim <- 1000
haiti_pop = 10911819
t_vacc_start <- dateToYears(as.Date('2019-01-12'))

time_frames = c(
  3,
  5,
  10)

scenarios = c(
  'S0',
  'S1',
  'S2',
  'S3',
  'S4',
  # 'S5',
  # 'S6',
  # 'S7',
  # 'S8',
  # 'S9',
  # 'S10',
  # 'S11',
  # 'S12',
  # 'S13',
  # 'S14',
  # 'S15',
  # 'S16',
  # 'S17',
  # 'S18',
  # 'S19',
  # 'S20',
  # 'S21',
  # 'S22',
  # 'S23',
  # 'S24',
  'S25'
  # 'S26',
  # 'S27',
  # 'S28',
  # 'S29',
  # 'S30',
  # 'S31',
  # 'S32',
  # 'S33',
  # 'S34',
  # 'S35',
  # 'S36'
)

thresh1 = haiti_pop / 10000
thresh2 = 1

pr_elim <- tribble(
  ~team, ~scenario, ~time_frame, ~pr_elim_thresh1, 
  ~pr_resurg_thresh1, ~pr_elim_thresh2, ~pr_resurg_thresh2,
  ~cum_inf_med, ~cum_inf_low, ~cum_inf_hi)

tte <-  tribble(
  ~team, ~scenario, ~t_elim_thresh1_med, ~t_elim_thresh1_low,
  ~t_elim_thresh1_hi, ~t_elim_thresh2_med, ~t_elim_thresh2_low, 
  ~t_elim_thresh2_hi, ~vacc_startdate)

ts <-  tribble(
  ~team, ~scenario, ~saturday_date,	~obs_cases_med, 
  ~obs_cases_low, ~obs_cases_hi, ~true_inf_med, ~true_inf_low,
  ~true_inf_hi, ~population_size, ~vacc_med, ~vacc_low, ~vacc_hi)

elimdate <-  tribble(
  ~team, ~scenario, ~elimination_date)


for (scenario in scenarios)
{
  print(scenario)
  load(file = sprintf("%sHaiti_OCV_Projection-allDep-%i-%s.rda", folder, nsim, scenario ))
  
  df <- projec %>% filter(sim != 'data', time > t_vacc_start) %>% select(time, sim, IncidenceAll)
  
  quantiles <- projec %>% select(time, sim, CasesAll, IncidenceAll, DosesAll) %>% as_tibble() %>% 
    mutate(isdata = sim == "data") %>%
    gather(variable, value, -time, -sim, -isdata) %>% 
    group_by(time, isdata, variable) %>% 
    summarise( q05 = quantile(value, 0.025, na.rm = T),
               mean = mean(value, na.rm = T),
               q50 = quantile(value, 0.5, na.rm = T),
               q95 = quantile(value, 0.975, na.rm = T)) %>%
    ungroup() %>%
    mutate(isdata = ifelse(isdata, "data", "simulation"),
           date = yearsToDateTime(time)) %>% 
    filter(date >= yearsToDate(sirb_cholera@t0)) %>%
    mutate(date = as.Date(round_date(date)))
  
  ts_inc <- quantiles %>% 
    mutate(saturday_date = as.Date(round_date(date))) %>% 
    filter(variable == "IncidenceAll", isdata == "simulation") %>% 
    select(saturday_date, q05, q50, q95) %>% 
    filter(dateToYears(saturday_date) > t_vacc_start) %>% 
    mutate(true_inf_med = q50) %>% 
    mutate(true_inf_low = q05) %>% 
    mutate(true_inf_hi = q95) %>%
    select(saturday_date, true_inf_med, true_inf_low, true_inf_hi)
  
  ts_obs <- quantiles %>% 
    mutate(saturday_date = as.Date(round_date(date))) %>% 
    filter(variable == "CasesAll", isdata == "simulation") %>% 
    select(saturday_date, q05, q50, q95) %>% 
    filter(dateToYears(saturday_date) > t_vacc_start) %>% 
    mutate(obs_cases_med = q50) %>% 
    mutate(obs_cases_low = q05) %>% 
    mutate(obs_cases_hi = q95) %>%
    select(saturday_date, obs_cases_med, obs_cases_low, obs_cases_hi)
  
  ts_vacc <- quantiles %>% 
    mutate(saturday_date = as.Date(round_date(date))) %>% 
    filter(variable == "DosesAll", isdata == "simulation") %>% 
    select(saturday_date, q05, q50, q95) %>% 
    filter(dateToYears(saturday_date) > t_vacc_start) %>% 
    mutate(vacc_med = q50) %>% 
    mutate(vacc_low = q05) %>% 
    mutate(vacc_hi = q95) %>%
    select(saturday_date, vacc_med, vacc_low, vacc_hi)
  
  dummy <- left_join(left_join(ts_inc, ts_obs), ts_vacc) %>% 
    mutate(population_size = haiti_pop) %>% 
    mutate(scenario = scenario) %>% 
    mutate(team = team) %>% 
    select(team, scenario, saturday_date,
           obs_cases_med, obs_cases_low, obs_cases_hi, 
           true_inf_med, true_inf_low, true_inf_hi,
           population_size, vacc_med, vacc_low, vacc_hi)
  
  ts <- rbind(ts, dummy)
  

  t_elim_thresh1 = list()
  t_elim_thresh2 = list()
  
  for (time_frame in time_frames)
  {
    print(time_frame)
    acc_elim_thresh1 = 0
    acc_elim_thresh2 = 0
    acc_resurg_thresh1 = 0
    acc_resurg_thresh2 = 0
    cum_inf = list()

    for (s in 1:nsim)
    {
      df_s <- df %>% filter(sim == s)
      roll_sum <- zoo::rollapply(df_s$IncidenceAll, 52, sum)
      
      if (sum(roll_sum[((time_frame-1)*52):(time_frame*52)] < thresh1, na.rm = T) >= 1 ){
        acc_elim_thresh1 = acc_elim_thresh1 + 1
        if (sum(df_s$IncidenceAll[(time_frame*52):length(roll_sum)] > thresh1/52.14, na.rm = T) >= 1 ){
          acc_resurg_thresh1 = acc_resurg_thresh1 +1
        }
      }
      if (sum(roll_sum[((time_frame-1)*52):(time_frame*52)] < thresh2, na.rm = T) >= 1 ){
        acc_elim_thresh2 = acc_elim_thresh2 + 1
        if (sum(df_s$IncidenceAll[(time_frame*52):length(roll_sum)] > thresh2/52.14, na.rm = T) >= 2 ){
          acc_resurg_thresh2 = acc_resurg_thresh2 +1  # TODO not length(roll_sum
        }
      }
      cum_inf <- c(cum_inf, df_s %>%
                     filter(time < t_vacc_start+time_frame) %>%
                     select(IncidenceAll) %>%
                     sum() )
      # Following code is time_frame independant:
      # Find if elimination
      if (time_frame == 3) {
        if (sum(roll_sum < thresh1, na.rm = T) >= 1 ){
          t_elim_thresh1 <- c(t_elim_thresh1, match(TRUE, roll_sum < thresh1))
        }
        if (sum(roll_sum < thresh2, na.rm = T) >= 1 ){
          t_elim_thresh2 <- c(t_elim_thresh2, match(TRUE, roll_sum < thresh2))
        }
        # Find elimination date
        if (sum(roll_sum < thresh2, na.rm = T) >= 1 ){
          ind <- match(TRUE, roll_sum < thresh2)
          elimination_date = toString(as.Date(yearsToDate(df_s$time[ind])))
          elimdate <- add_row(elimdate, team, scenario, elimination_date)
        } else {
          elimination_date = NA
          elimdate <- add_row(elimdate, team, scenario, elimination_date)
        }
      }
    }
    
    
    pr_elim_thresh1 =   acc_elim_thresh1 / nsim
    pr_elim_thresh2 =   acc_elim_thresh2 / nsim
    pr_resurg_thresh1 =   acc_resurg_thresh1 / acc_elim_thresh1
    pr_resurg_thresh2 =   acc_resurg_thresh2 / acc_elim_thresh2
    
    cum_inf <- unlist(cum_inf)
    
    cum_inf_med = quantile(cum_inf, 0.5, na.rm = T)
    cum_inf_low = quantile(cum_inf, 0.025, na.rm = T)
    cum_inf_hi = quantile(cum_inf, 0.975, na.rm = T)
    
    pr_elim <- add_row(pr_elim, team, scenario, time_frame, pr_elim_thresh1, 
                       pr_resurg_thresh1, pr_elim_thresh2, pr_resurg_thresh2,
                       cum_inf_med, cum_inf_low, cum_inf_hi)
  }
  
  t_elim_thresh1 <- unlist(t_elim_thresh1)
  t_elim_thresh2 <- unlist(t_elim_thresh2)
  t_elim_thresh1_med = quantile(t_elim_thresh1, 0.5, na.rm = T)
  t_elim_thresh1_low = quantile(t_elim_thresh1, 0.025, na.rm = T)
  t_elim_thresh1_hi  = quantile(t_elim_thresh1, 0.975, na.rm = T)
  t_elim_thresh2_med = quantile(t_elim_thresh2, 0.5, na.rm = T)
  t_elim_thresh2_low = quantile(t_elim_thresh2, 0.025, na.rm = T)
  t_elim_thresh2_hi  = quantile(t_elim_thresh2, 0.975, na.rm = T)
  vacc_startdate = '2019-01-12'
  
  tte <-  add_row(tte,
    team, scenario, t_elim_thresh1_med, t_elim_thresh1_low,
    t_elim_thresh1_hi, t_elim_thresh2_med, t_elim_thresh2_low, 
    t_elim_thresh2_hi, vacc_startdate)
  
}

write.csv(pr_elim, file = paste0(folder, "pr_elim_", team, '_', Sys.Date(), '.csv'), row.names=FALSE)
write.csv(ts, file = paste0(folder, "ts_", team, '_', Sys.Date(), '.csv'), row.names=FALSE)
write.csv(tte, file = paste0(folder, "tte_", team, '_', Sys.Date(), '.csv'), row.names=FALSE)
write.csv(elimdate, file = paste0(folder, "elimdate_", team, '_', Sys.Date(), '.csv'), row.names=FALSE)

