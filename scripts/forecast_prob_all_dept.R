library(tidyverse)
library(magrittr)
library(ggthemes)
library(GGally)
library(foreach)
library(itertools)
library(pomp)

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
scenario <- 'S25'
haiti_pop = 10911819
t_vacc_start <- dateToYears(as.Date('2019-01-12'))

time_frames = c(
  3,
  5,
  10)

scenarios = c(
 # 'S0',
  # 'S1',
  #'S2',
  #'S3',
  #'S4',
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

ts <- tribble(
  ~team, ~scenario, ~saturday_date,
  ~obs_cases_med, ~obs_cases_low, ~obs_cases_hi, 
  ~true_inf_med, ~true_inf_low, ~true_inf_hi,
  ~population_size, ~vacc_med, ~vacc_low, ~vacc_hi)

for (scenario in scenarios)
{
  load(file = sprintf("%sHaiti_OCV_Projection-allDep-%i-%s.rda", folder, nsim, scenario ))
  
  
  for (time_frame in time_frames)
  {
    acc_elim_thresh1 = 0
    acc_elim_thresh2 = 0
    acc_resurg_thresh1 = 0
    acc_resurg_thresh2 = 0
    
    for (s in 1:nsim) {
      df_s <- df %>% filter(sim == s)
      if (sum(
        df_s %>% filter(time > t_vacc_start + 3 - 6 / 12)
        %>% filter(time < t_vacc_start + 3 + 6 / 12)
        %>% select('sum')
      ) < thresh1) {
        acc_elim_thresh1 = acc_elim_thresh1 + 1
      }
      pr_elim_thresh1 =   acc_elim_thresh1 / nsim
      pr_elim_thresh2 =   acc_elim_thresh2 / nsim
      pr_resurg_thresh1 =   acc_resurg_thresh1 / nsim
      pr_resurg_thresh2 =   acc_resurg_thresh2 / nsim
      
    pr_elim <- add_row(pr_elim, team, scenario, time_frame, pr_elim_thresh1, 
                       pr_resurg_thresh1, pr_elim_thresh2, pr_resurg_thresh2,
                       cum_inf_med, cum_inf_low, cum_inf_hi)
  }
}

  
    if (dp == 'Artibonite') {
      df <-
        projec %>% filter(sim != 'data', time > t_vacc_start) %>% select(time, sim, I)
      colnames(df)[names(df) == "I"] <- dp
      
    } else{
      dat <-
        projec %>% filter(sim != 'data', time > 2019) %>% select(time, sim, I)
      colnames(dat)[names(dat) == "I"] <- dp
      
      df <- merge(df, dat, by = c('time', 'sim'))
    }
  }
  
  df$sum <- rowSums(df %>% select(-sim,-time))
  
  df %<>% select(sim, time, sum)

  

    if (sum(
      df_s %>% filter(time > t_vacc_start + 3 - 6 / 12)
      %>% filter(time < t_vacc_start + 3 + 6 / 12)
      %>% select('sum')
    ) < thresh2) {
      acc_3y_thresh2 = acc_3y_thresh2 + 1
    }
    
    if (sum(
      df_s %>% filter(time > t_vacc_start + 5 - 6 / 12)
      %>% filter(time < t_vacc_start + 5 + 6 / 12)
      %>% select('sum')
    ) < thresh1) {
      acc_5y_thresh1 = acc_5y_thresh1 + 1
    }
    if (sum(
      df_s %>% filter(time > t_vacc_start + 5 - 6 / 12)
      %>% filter(time < t_vacc_start + 5 + 6 / 12)
      %>% select('sum')
    ) < thresh2) {
      acc_5y_thresh2 = acc_5y_thresh2 + 1
    }
    
    
    if (sum(
      df_s %>% filter(time > t_vacc_start + 10 - 6 / 12)
      %>% filter(time < t_vacc_start + 10 + 6 / 12)
      %>% select('sum')
    ) < thresh1) {
      acc_10y_thresh1 = acc_10y_thresh1 + 1
    }
    if (sum(
      df_s %>% filter(time > t_vacc_start + 10 - 6 / 12)
      %>% filter(time < t_vacc_start + 10 + 6 / 12)
      %>% select('sum')
    ) < thresh2) {
      acc_10y_thresh2 = acc_10y_thresh2 + 1
    }
  }
  
  print(paste0('3 years T1 : ', acc_3y_thresh1 / nsim * 100))
  print(paste0('3 years T2 : ', acc_3y_thresh2 / nsim * 100))
  
  print(paste0('5 years T1 : ', acc_5y_thresh1 / nsim * 100))
  print(paste0('5 years T2 : ', acc_5y_thresh2 / nsim * 100))
  
  print(paste0('10 years T1 : ', acc_10y_thresh1 / nsim * 100))
  print(paste0('10 years T2 : ', acc_10y_thresh2 / nsim * 100))
  
  probability_data[scenario,] = c(acc_3y_thresh1 / nsim * 100,
                                  acc_3y_thresh2 / nsim * 100,
                                  acc_5y_thresh1 / nsim * 100,
                                  acc_5y_thresh2 / nsim * 100,
                                  acc_10y_thresh1 / nsim * 100,
                                  acc_10y_thresh2 / nsim * 100)
  
}

save(probability_data,file="probability_data.Rda")

