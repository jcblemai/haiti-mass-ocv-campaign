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
coef(sirb_cholera)["Rtot_0"] <-25

p <- simulate(sirb_cholera, nsim = 10, as.data.frame = T) %>% 
  gather(variable, value, -time, -rain, -sim)  %>% 
  ggplot(aes(x = time, y = value, color = sim)) + 
  #  geom_line(aes(y = cases), color = datacol, lwd = 0.2) 
  geom_line() + 
  facet_wrap(~variable, scales = "free_y") 


datacol <- "#ED0000"
print(p)

# several simulation --------------------------------------------------------------


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

yr_thresh <- 2015
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
  
