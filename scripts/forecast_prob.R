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

# TODO comment to run in python

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
nsim <- 500
scenario <- 'S25'
haiti_pop = 10911819
t_vacc_start <- dateToYears(as.Date('2019-01-12'))

scenarios = c(
  'S1',
  'S2',
  'S3',
  'S4',
  'S5',
  'S6',
  'S7',
  'S8',
  'S9',
  'S10',
  'S11',
  'S12',
  'S13',
  'S14',
  'S15',
  'S16',
  'S17',
  'S18',
  'S19',
  'S20',
  'S21',
  'S22',
  'S23',
  'S24',
  'S25',
  'S26',
  'S27',
  'S28',
  'S29',
  'S30',
  'S31',
  'S32',
  'S33',
  'S34',
  'S35',
  'S36'
)


probability_data = data.frame(scenarios, 'y3T1' = 0,  'y3T2' = 0, 
                                           'y5T1' = 0,  'y5T2' = 0,
                                           'y10T1' = 0, 'y10T2' = 0, row.names = 1)


probability_data


for (scenario in scenarios)
{
  for (dp in departements) {
    print(dp)
    load(
      file = sprintf(
        "output/Simulations/Haiti_OCV_Projection-%s-%i-%s.rda",
        dp,
        nsim,
        scenario
      )
    )
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
  
  thresh1 = haiti_pop / 10000
  thresh2 = 1
  
  acc_3y_thresh1 = 0
  acc_3y_thresh2 = 0
  
  acc_5y_thresh1 = 0
  acc_5y_thresh2 = 0
  
  acc_10y_thresh1 = 0
  acc_10y_thresh2 = 0
  
  
  
  
  for (s in 1:nsim) {
    df_s <- df %>% filter(sim == s)
    if (sum(
      df_s %>% filter(time > t_vacc_start + 3 - 6 / 12)
      %>% filter(time < t_vacc_start + 3 + 6 / 12)
      %>% select('sum')
    ) < thresh1) {
      acc_3y_thresh1 = acc_3y_thresh1 + 1
    }
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

