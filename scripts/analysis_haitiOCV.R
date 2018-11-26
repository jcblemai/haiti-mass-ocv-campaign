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
}
departement <- args[1]


# Pair plots ---------------------------------------------------------------

# Stochastic model (POMP)
liks_stoch <- read_csv("results/Haiti_OCVparam_logliks-10-l1.csv")
liks <- liks_stoch

#colnames(liks_stoch) <- str_replace_all(colnames(liks_stoch), "_", "")

doplots <- T

if(doplots) {
  
  plotPairs <- function(data, variables, filename, width = 12, height = 8) {
    p <- ggpairs(data %>% 
                   select(loglik, one_of(variables)) %>% 
                   keep(~sd(.) > 1e-4) %>% 
                   map_df(~ map_dbl(., ~ifelse(. == 0, NA, .))), #%>% bind_cols(data[,"model"]),
                 # filter(loglik > -400),
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
      filename = str_c("results/figures/all-logliks.png"),
      width = 10, height = 8)
  
  # plot all
  plotPairs(liks_stoch %>% 
              arrange(desc(loglik)) %>% 
              slice(1:20) %>%
              ungroup,
            c("sigma", "betaB", "mu_B" , "thetaA", "thetaI", "rhoA", "rhoI", "r", "lambda"),
            "results/figures/stoch_posteriors.png",
            width = 12,
            height = 12)
  
  # plot pairs by type of parameteres
  plotPairs(liks_stoch, c("epsilon", "k"), "results/figures/liks_measurement_model.png", width = 10, height = 6.5)
  plotPairs(liks_stoch, str_c(c("S", "A", "I", "R", "B"), "_0"), "results/figures/liks_initial_conditions.png", width = 10, height = 5)
  plotPairs(liks_stoch, c("sigma", "betaB", "mu_B", "thetaA", "thetaI", "rhoA", "rhoI",  "std_W"), "results/figures/liks_sirb_processes.png")
  plotPairs(liks_stoch, c("sigma", "betaB", "mu_B", "thetaA", "thetaI", "lambda", "r"), "results/figures/liks_rainfall_effect.png")

}
# Likelihood comparison ---------------------------------------------------
load(paste0("sirb_cholera_pomped_", departement, ".rda"))


# get best likelihood
best_liks <- liks %>% 
arrange(desc(loglik)) %>% 
slice(1)
  

# Get estimated parameter names
param_est_names_stoch <- colnames(liks_stoch)  %>% 
  keep(~ !(. %in% c("loglik", "loglik.se", "type", "k")) & sd(liks_stoch[[.]]) > 1e-3)


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
best_param <- best_liks %>% 
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
sim_stochastic <- foreach(r = iter(best_param, by = "row"), 
                          .combine = rbind) %dopar% {
                            
                          simulatePOMP(sirb_cholera, unlist(r[names(coef(sirb_cholera))]), nsim = 10)
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
  scale_x_date(date_labels = "%b-%y", expand = c(0,0), limits = as.Date(c("2015-06-01", "2018-01-01"))) +
  scale_y_continuous(expand = c(0,0))+ 
  labs(y = "daily cholera cases", x = "date") +
  theme(panel.grid.major = element_line(color = "lightgray"),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        strip.text = element_blank(),
        axis.title = element_text())

p.sim
ggsave(p.sim, filename = "results/figures/simulations_comparison.png", width = 9, height = 10, dpi = 300)


# Parameter profiles ------------------------------------------------------
# function to compute confidence intervales from Ionides et al., 2017
mcap <- function(lp,parameter,confidence=0.95,lambda=0.75,Ngrid=1000){
  smooth_fit <- loess(lp ~ parameter, span=lambda)
  parameter_grid <- seq(min(parameter), max(parameter), length.out = Ngrid)
  smoothed_loglik <- predict(smooth_fit,newdata=parameter_grid)
  smooth_arg_max <- parameter_grid[which.max(smoothed_loglik)]
  smooth_ll_max <- max(smoothed_loglik)
  
  dist <- abs(parameter-smooth_arg_max)
  included <- dist < sort(dist)[trunc(lambda*length(dist))]
  maxdist <- max(dist[included])
  weight <- rep(0,length(parameter))
  weight[included] <- (1-(dist[included]/maxdist)^3)^3
  quadratic_fit <- lm(lp ~ a + b, weight=weight,
                      data = data.frame(lp=lp,b=parameter,a=-parameter^2)
  )
  b <- unname(coef(quadratic_fit)["b"] )
  a <- unname(coef(quadratic_fit)["a"] )
  m <- vcov(quadratic_fit)
  var_b <- m["b","b"]
  var_a <- m["a","a"]
  cov_ab <- m["a","b"]
  se_mc_squared <- (1 / (4 * a^2)) * (var_b - (2 * b/a) * cov_ab + (b^2 / a^2) * var_a)
  se_stat_squared <- 1/(2*a)
  se_total_squared <- se_mc_squared + se_stat_squared
  delta <- qchisq(confidence,df=1) * ( a * se_mc_squared + 0.5)
  loglik_diff <- max(smoothed_loglik) - smoothed_loglik
  ci <- range(parameter_grid[loglik_diff < delta])
  
  list(lp=lp,parameter=parameter,confidence=confidence,
       quadratic_fit=quadratic_fit, quadratic_max=b/(2*a),
       smooth_fit=smooth_fit,
       fit=data.frame(
         parameter=parameter_grid,
         smoothed=smoothed_loglik,
         quadratic=predict(quadratic_fit, list(b = parameter_grid, a = -parameter_grid^2))
       ),
       llmax = smooth_ll_max,
       mle=smooth_arg_max, ci=ci, delta=delta,
       se_stat=sqrt(se_stat_squared), se_mc=sqrt(se_mc_squared), se=sqrt(se_total_squared)
  )
} 


profile_liks <- foreach(param = c("alpha_E", "lambda_E"), .combine = rbind) %do% {
  load(str_c("results/choleraJuba-SIR_B_E-", param ,"-2-profiles_3.rda"))
  liks <- lik_mf
  load(str_c("results/choleraJuba-SIR_B_E-", param ,"-2-profiles_4.rda"))
  liks <- rbind(liks, lik_mf)
  
  liks %>%  
    select(loglik, one_of(param)) %>% 
    gather(variable, value, -loglik)
} %>% 
  filter(loglik > -330)# %>%
# group_by(variable, value) %>%
# arrange(desc(loglik)) %>%
# slice(1) %>%
# ungroup 


removeOutilers <- function(df){
  ll.sd <- sd(df[["loglik"]])
  ll.m <- mean(df[["loglik"]])
  
  df[abs(df[["loglik"]] - ll.m) < (1 * ll.sd) , ]
}

# profile_liks %<>% 
#   group_by(variable) %>% 
#   mutate(cut = cut(value, breaks = 15)) %>% 
#   ungroup %>% 
#   nest(-variable, -cut) %>% 
#   mutate(data2 = map(data, ~removeOutilers(.))) %>% 
#   select(-cut,-data) %>% 
#   unnest

# compute likelihood profiles
param_mcap <- profile_liks %>% 
  nest(-variable) %>% 
  mutate(mcap = map(data, ~mcap(.$loglik, .$value, lambda = 1))) 

param_cis <- param_mcap %>% 
  mutate(ci = map(mcap, "ci")) %>% 
  select(-data, -mcap) %>% 
  unnest() %>% 
  group_by(variable) %>% 
  mutate(bound = ifelse(ci == min(ci), "lower","upper")) %>% 
  ungroup %>% 
  spread(bound, ci)

param_mle <- param_mcap %>% 
  mutate(vec = map(mcap, ~ c(.$llmax, .$delta))) %>% 
  select(-data, -mcap) %>% 
  unnest() %>% 
  group_by(variable) %>% 
  mutate(bound = ifelse(vec < 0, "llmax", "delta")) %>% 
  ungroup %>% 
  spread(bound, vec)

param_quad <- param_mcap %>% 
  mutate(fit = map(mcap, "fit")) %>% 
  select(-data, -mcap) %>% 
  unnest() %>% 
  gather(fun, value, -variable, -parameter) %>% 
  mutate(fun = factor(fun, levels = c("smoothed", "quadratic"))) 


# plot
p.profile <- profile_liks %>% 
  ggplot(aes(x = value, y = loglik)) +
  geom_line(data = param_quad , aes(x = parameter, y = value, color = fun, linetype = fun)) +
  guides(color = "none", linetype = "none") +
  geom_vline(data = param_cis, aes(xintercept = lower), color = "red") +
  geom_vline(data = param_cis, aes(xintercept = upper), color = "red") +
  geom_hline(data = param_mle, aes(yintercept = llmax - delta), color = "red") +
  geom_point() +
  scale_color_manual(values = c("red", "blue")) +
  theme_bw() +
  facet_wrap(~variable, scales = "free_x")

ggsave(p.profile, filename = "results/figures/profiles_alpha_lambda_E.png", width = 10, height = 5, dpi = 300)