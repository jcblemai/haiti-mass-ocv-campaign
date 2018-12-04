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

# One simulation  --------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  # default departement
  args[1] = "Artibonite"
}
departement <- args[1]

load(paste0("sirb_cholera_pomped_", departement, ".rda")
     )


coef(sirb_cholera)["mu_B"] <- 365/5
coef(sirb_cholera)["betaB"] <- 0.1
coef(sirb_cholera)["Rtot_0"] <-50

p <- simulate(sirb_cholera, nsim = 10, as.data.frame = T) %>% 
  gather(variable, value, -time, -rain, -sim)  %>% 
  ggplot(aes(x = time, y = value, color = sim)) + 
  #  geom_line(aes(y = cases), color = datacol, lwd = 0.2) 
  geom_line() + 
  facet_wrap(~variable, scales = "free_y") 


datacol <- "#ED0000"
print(p)

# several simulation --------------------------------------------------------------


doMC::registerDoMC(8)

# function to simulate from a given set of paramters
simulatePOMP <- function(pmp, params, nsim, seed = 199919L) {
  pmpsim <- pmp
  pomp::coef(pmpsim) <- params
  
  pomp::simulate(pmpsim, nsim = nsim, as.data.frame = T , include.data = TRUE, seed = seed) %>% 
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
    filter(date >= yearsToDate(pmpsim@t0))
}

# run simulations for each model
nsim = 10
sim_stochastic <- foreach(r = iter(best_param, by = "row"), #over models, useless
                          .combine = rbind) %do% {
                            
                            simulatePOMP(sirb_cholera, unlist(r[names(coef(sirb_cholera))]), nsim = nsim)
                          }

# tidy tibble for merger
sim_stochastic_quantiles <- sim_stochastic %>% 
  mutate(date = as.Date(round_date(date)), type = "stochastic") %>% 
  filter(variable == "cases", isdata == "simulation") %>% 
  select(-isdata, -time, -variable) %>% 
  bind_cols(sim_stochastic %>% 
              filter(variable == "cases", isdata == "data") %>% 
              select(mean) %>% 
              rename(cases = mean))

# PLOT BOTH
simcol <- "#175CD6"
datacol <- "#ED0000"

p.sim <- ggplot(data = sim_stochastic_quantiles,
                aes(x = date))+
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.1, color = simcol, fill = simcol) +
  geom_line(aes(y = q50), color = simcol) +
  geom_line(aes(y = mean), linetype = 2, color = simcol) +
  geom_line(aes(y = cases), color = datacol, lwd = 0.2) +
  geom_point(aes(y = cases), color = datacol, size = 0.8) +
  #geom_text(data = psim_labels, aes (y = y, label = label), size = 7) +
  #facet_grid(model~type) +
  scale_x_date(date_labels = "%b-%y", expand = c(0,0), limits = as.Date(c("2014-03-01", "2018-07-14"))) +
  scale_y_continuous(expand = c(0,0))+ 
  labs(y = "daily cholera cases", x = "date") +
  theme(panel.grid.major = element_line(color = "lightgray"),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        strip.text = element_blank(),
        axis.title = element_text())

p.sim


# Data plots --------------------------------------------------------------

cases <- read_csv("haiti-data/fromAzman/cases_corrected.csv")  %>% 
  gather(dep, cases, -date) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date))


rain <- read_csv("haiti-data/fromAzman/rainfall.csv")  %>% 
  gather(dep, rain, -date) %>% 
  group_by(dep) %>% 
  mutate(max_rain = max(rain), rain_std = rain/max_rain) %>%
  ungroup() %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date))


library(tidyverse)

yr_thresh <- 2012
cases %>% 
  mutate(yr = lubridate::year(date)) %>%
  filter(yr>=yr_thresh) %>% 
  mutate(incal = case_when(date > as.Date("2014-03-01") ~ T,
                           T ~ F)) %>% 
  ggplot(aes(x = date, y = cases)) +
  geom_bar(data = filter(rain, lubridate::year(date)>=yr_thresh), aes(y = rain * 10), fill = "blue", stat = "identity", alpha = 0.4) +
  geom_line(aes(color = incal)) +
  scale_color_manual(values = c("black", "red")) +
  # geom_vline(aes(xintercept = as.Date("2017-07-08"))) + 
  # geom_vline(aes(xintercept = as.Date("2018-07-14"))) + 
  coord_cartesian(ylim = c(0,500)) +
  facet_wrap(~dep)


# Data plots togetehr  --------------------------------------------------------------

cases <- read_csv("haiti-data/fromAzman/cases_corrected.csv")  %>% 
  gather(dep, cases, -date) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date))


rain <- read_csv("haiti-data/fromAzman/rainfall.csv")  %>% 
  gather(dep, rain, -date) %>% 
  group_by(dep) %>% 
  mutate(max_rain = max(rain), rain_std = rain/max_rain) %>%
  ungroup() %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         time = dateToYears(date))


library(tidyverse)

yr_thresh <- 2011
cases %>% 
  mutate(yr = lubridate::year(date)) %>%
  filter(yr>=yr_thresh) %>% 
  mutate(incal = case_when(date > as.Date("2016-07-09") ~ T,
                           T ~ F)) %>% 
  ggplot(aes(x = date, y = cases)) +
  geom_bar(data = filter(rain, lubridate::year(date)>=yr_thresh), aes(y = rain * 10), fill = "blue", stat = "identity", alpha = 0.4) +
  geom_line(aes(color = incal)) +
  scale_color_manual(values = c("black", "red")) +
  # geom_vline(aes(xintercept = as.Date("2017-07-08"))) + 
  # geom_vline(aes(xintercept = as.Date("2018-07-14"))) + 
  coord_cartesian(ylim = c(0,1500)) +
  facet_wrap(~dep)



cases %>% 
  mutate(yr = lubridate::year(date)) %>%
  group_by(yr) %>% 
  mutate(julian = julian(date, min(date))) %>% 
  filter(yr>=yr_thresh) %>%
  filter(dep == "Artibonite") %>% 
  # mutate(incal = case_when(date > as.Date("2016-07-09") ~ T,
  # T ~ F)) %>% 
  ggplot(aes(x = julian, y = cases)) +
  geom_bar(data = filter(rain, dep == "Artibonite", year(date)>=yr_thresh) %>% 
             mutate(yr = year(date)) %>% 
             group_by(yr) %>% 
             mutate(julian = julian(date, min(date))) ,
           aes(y = rain * 10), stat = "identity", alpha = 0.4) +
  geom_line() +
  facet_grid(yr~.)
  
