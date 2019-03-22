# Title: Mif run results
# Description: Explore results of mif runs
# Date: Fri Jul 13 10:55:26 2018
# Author: javier.perezsaez@epfl.ch



# Preamble ---------------------------------------------------------------
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

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  # default departement
  args[1] = "Artibonite"
  args[2] = 4
} else if (length(args)==1) {
  args[2] = 4
}


departement <- args[1]
run_level <- as.integer(args[2])

departement <- 'Artibonite'
run_level <- 4
output_dir <- "output_15-05-long/"


# Pair plots ---------------------------------------------------------------
# One line per initial condition
liks_stoch <- read_csv(sprintf("%s%s/Haiti_OCV-%s-param_logliks-10-l%i.csv", output_dir, departement, departement, run_level)) 
liks <- liks_stoch

#colnames(liks_stoch) <- str_replace_all(colnames(liks_stoch), "_", "")

doplots <- T

if(doplots) {
  
  plotPairs <- function(data, variables, filename, width = 12, height = 8) {
    p <- ggpairs(data %>% 
                   select(loglik, one_of(variables)) %>% 
                   keep(~sd(.) > 1e-4) %>% 
                   map_df(~ map_dbl(., ~ifelse(. == 0, NA, .))), #%>% 
#                  filter(loglik > -750),
                 aes(alpha = I(0.4)),
                 upper = list(continuous = "points", combo = "box_no_facet", discrete = "facetbar", na = "na"),
                 lower = list(continuous = "points", combo = "facethist", discrete = "facetbar", na = "na") ,
                 legend = NULL
    ) +
      theme_few() +
      # scale_fill_few()+
      theme(
        # panel.grid.major = element_line(color = "grey"),
        # axis.line = element_line(color = "black"),
        panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = .3)
      ) 
    
    p
    ggsave(p, filename = filename, width = width, height = width)
  }
  
  # create pairplots by model
  ggsave(ggpairs(liks_stoch %>% select(-loglik.se, -contains("0")) %>% keep(~sd(.) > 1e-4),
              upper = list(continuous = "points", combo = "box_no_facet", discrete = "facetbar", na = "na"),
              lower = list(continuous = "cor", combo = "facethist", discrete = "facetbar", na = "na"),
              title = sprintf("Haiti OCV model with best loglik: %f",  max(liks_stoch$loglik))
      ),
      filename = str_c(output_dir, departement, "/all-logliks.png"),
      width = 10, height = 8)
  
  # plot all
  plotPairs(liks_stoch %>% 
              arrange(desc(loglik)) %>% 
              slice(1:20) %>%
              ungroup,
            c("sigma", "betaB", "mu_B" , "XthetaA", "thetaI", "rhoA", "XrhoI", "r", "lambda"),
            str_c(output_dir, departement, "/stoch_posteriors.png"),
            width = 12,
            height = 12)
  
  # plot pairs by type of parameteres
  plotPairs(liks_stoch, c("epsilon", "k"),  str_c(output_dir, departement, "/liks_measurement_model.png"), width = 10, height = 6.5)
  plotPairs(liks_stoch, c("Rtot_0"),  str_c(output_dir, departement, "/liks_initial_conditions.png"), width = 10, height = 5)
  plotPairs(liks_stoch, c("sigma", "betaB", "mu_B", "XthetaA", "thetaI", "rhoA", "XrhoI",  "std_W"),  str_c(output_dir, departement, "/liks_sirb_processes.png"))
  plotPairs(liks_stoch, c("sigma", "betaB", "mu_B", "XthetaA", "thetaI", "lambda", "r", "Rtot_0"),  str_c(output_dir, departement, "/liks_rainfall_effect.png"))

}
# Likelihood comparison ---------------------------------------------------
load(paste0(output_dir, departement, "/sirb_cholera_pomped_", departement, ".rda"))

# Compare outputs ---------------------------------------------------------

doMC::registerDoMC(8)

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


# get MLE paramter sets 
best_param <- liks %>% 
arrange(desc(loglik)) %>% 
slice(1)  %>% 
  arrange(desc(loglik)) %>% 
  ungroup %>% 
  left_join(liks_stoch)

# How to change covariate if necessary: Javier recrated a pmp object with a covariate he can change
# re-create pomp object for simulations
# sirb_cholera_maxrain <- pomp(sirb_cholera,
#                              covar = rain %>%
#                                filter(time >= min(sirb_cholera@tcovar) & time <= max(sirb_cholera@tcovar)) %>% 
#                                select(time, rain_std) %>%
#                                rename(rain = rain_std),
#                              tcovar = "time")

# run simulations for each model
nsim = 1000
sim_stochastic <- simulatePOMP(sirb_cholera, unlist(best_param[names(coef(sirb_cholera))]), nsim = nsim)

# tidy tibble for merger
sim_stochastic_quantiles <- sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "cases", isdata == "simulation") %>% 
  select(-isdata, -time, -variable) %>% 
  bind_cols(sim_stochastic %>%                                         #Add a column with the data
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
ggsave(p.sim, filename = str_c(output_dir, departement, "/simulations_comparison.png"), width = 9, height = 10, dpi = 300)

sim_stochastic_quantiles_all <- sim_stochastic %>% 
  mutate(date = as.Date(round_date(date))) %>% 
  filter(variable == "cases" |  isdata == "simulation") %>% 
  select(-isdata) 

p.all <- ggplot(data = sim_stochastic_quantiles_all,
                aes(x = date))+
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.1, color = simcol, fill = simcol) +
  geom_line(aes(y = q50), color = simcol) +
  geom_line(aes(y = mean), linetype = 2, color = simcol) +
    facet_wrap(~variable, scales = "free_y") 
  scale_x_date(date_labels = "%b-%y", expand = c(0,0), limits = as.Date(c("2014-03-01", "2018-07-14")))


p.all
ggsave(p.all, filename = str_c(output_dir, departement, "/simulations_all.png"), width = 9, height = 10, dpi = 300)
