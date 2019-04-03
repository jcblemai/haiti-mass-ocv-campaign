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

# Analysis Calibration  --------------------------------------------------------------
output_dir <- "output_15-28-cdc/"
departement <- "Artibonite"
run_level <- 3

load(paste0(output_dir,departement, '/HaitiOCV-',run_level,'-',departement, '-mif_runs.rda' ))
plot(mf)

# Analysis Projection -----------------------------------------------------------
load('output/Simulations/Haiti_OCV_Projection-Artibonite-500-S0.rda')
load(paste0(output_dir, departement, "/sirb_cholera_pomped_", departement, ".rda"))

d <- projec %>%
  gather(variable, value, -time, -rain, -sim)  %>% 
  ggplot(aes(x = time, y = value, color = sim)) + 
  #  geom_line(aes(y = cases), color = datacol, lwd = 0.2) 
  geom_line() + 
  facet_wrap(~variable, scales = "free_y") +  
  theme(legend.position = "none") # Legend take too much space

datacol <- "#ED0000"
print(d)


# Analysis Projection -----------------------------------------------------------
departement <- 'Artibonite'
nsim <- 500
scenario <- 'S0'
load(paste0(output_dir, departement, "/sirb_cholera_pomped_", departement, ".rda"))
load(sprintf("output/Simulations/Haiti_OCV_Projection-%s-%i-%s.rda", departement, nsim, scenario))
input_parameters <- yaml::read_yaml("haiti-data/input_parameters.yaml")

# Start and end dates of epidemic
t_start <- dateToYears(as.Date(input_parameters$t_start))
t_end <- dateToYears(as.Date(input_parameters$t_end))
t_forecast <- dateToYears(as.Date("2029-12-21"))


quant <- projec %>%
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


simcol <- "#175CD6"
datacol <- "#ED0000"

p.sim <- ggplot(data = quant,
                 aes(x = date))+
   geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.1, color = simcol, fill = simcol) +
   geom_line(aes(y = q50), color = simcol) + 
  facet_wrap(~variable, scales = "free_y") +  
  theme(legend.position = "none") + # Legend take too much space
  #theme(panel.grid.major = element_line(color = "lightgray"),
  #      panel.background = element_rect(fill = "white"),
  #      axis.line = element_line(color = "black"),
  #      strip.text = element_blank(),
  #      axis.title = element_text()) +
#scale_y_continuous(expand = c(0,0))+
   labs(y = "daily cholera cases", x = "date") #+
  # theme(panel.grid.major = element_line(color = "lightgray"),
  #       panel.background = element_rect(fill = "white"),
  #       axis.line = element_line(color = "black"),
  #       strip.text = element_blank(),
  #         axis.title = element_text())
 
 print(p.sim)
# One simulation  --------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  # default departement
  args[1] = "Artibonite"
}
departement <- args[1]
output_dir <- "output/"
load(paste0(output_dir, departement, "/sirb_cholera_pomped_", departement, ".rda")
     )


coef(sirb_cholera)["mu_B"] <- 365/5
coef(sirb_cholera)["betaB"] <- 0.1
coef(sirb_cholera)["Rtot_0"] <-25

# Initialize the parameters to estimate (just initial guesses)
# param_est <- set_names(seq_along(param_est_names) * 0, param_est_names)
# param_est["sigma"] <- .2
# param_est["rhoA"] <- 1/(365*3)
# param_est["XrhoI"] <- 0.5
# param_est["betaB"] <- .1
# param_est["mu_B"] <-  365/5
# param_est["XthetaA"] <- 0.5
# param_est["thetaI"] <- .01
# param_est["lambda"] <- 0
# param_est["lambdaR"] <- 10
# param_est["r"] <- 1
# param_est["std_W"] <- .001
# param_est["epsilon"] <- .5
# param_est["k"] <- 10001
# param_est["Rtot_0"] <- 0.35
# param_est["foi_add"] <- 0.001
# param_est["gammaA"] <- 1/5
# m_est["gammaI"] <- 1/5

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

yr_thresh <- 2010
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
 # coord_cartesian(ylim = c(0,500)) +
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
  
