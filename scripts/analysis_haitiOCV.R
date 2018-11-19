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
# load fits ---------------------------------------------------------------

# Stochastic model (POMP)
liks_stoch <- read_csv("results/Haiti_OCVparam_logliks-10-l1.csv")

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
load("data/sirb_cholera_pomped.rda")

# function to compute nesting among models
isNested <- function(param_vec1, param_vec2) {
  isnested <-  sum(names(param_vec1)[param_vec1!=0] %>% 
                     map_dbl(~ sum(!(. %in% names(param_vec2)[param_vec2!=0])))) == 0
  
  dof_diff <-  sum(names(param_vec2)[param_vec2!=0] %>% 
                     map_dbl(~ sum(!(. %in% names(param_vec1)[param_vec1!=0]))))
  
  return(c(isnested = isnested, dof_diff = dof_diff, nparam_from = sum(param_vec1!=0), nparam_to = sum(param_vec2!=0)))
}


computeRatioTests <- function(liks, param_est_names, model.type)  {
  
  # get best likelihoods for each model
  best_liks <- liks %>% 
    group_by(model) %>% 
    arrange(desc(loglik)) %>% 
    slice(1)
  
  # determine model nestings based on paramter values
  model_nestings <- foreach(mod = iter(best_liks, by = "row"),
                            .combine = rbind) %do% {
                              foreach(mod2 = iter(best_liks %>%
                                                    filter(model != mod$model), by = "row"),
                                      .combine = rbind) %do% {
                                        data.frame(model_ref = mod$model,
                                                   model_comp = mod2$model,
                                                   as.list(isNested(mod[,param_est_names], mod2[,param_est_names]))) %>%
                                          mutate(deltaLL = 2* (mod2$loglik - mod$loglik)) %>% 
                                          as_tibble()
                                      }
                            } # %>
  model_nestings$LLtest = apply(model_nestings[, -c(1:2)] %>% as.matrix, 1, function(x) ifelse(x["deltaLL"] > 0, dchisq(x["deltaLL"], x["dof_diff"]), NA))
  model_nestings$LLtest2 = apply(model_nestings[, -c(1:2)] %>% as.matrix, 1, function(x) ifelse(x["deltaLL"] > 0, x["deltaLL"] > qchisq(0.95, x["dof_diff"]), NA))
  
  model_nestings$type = model.type
  return(model_nestings)
}

# Get estimated parameter names
param_est_names_stoch <- colnames(liks_stoch)  %>% 
  keep(~ !(. %in% c("loglik", "loglik.se", "model", "type", "k")) & sd(liks_stoch[[.]]) > 1e-3)

param_est_names_deter <- colnames(liks_deter) %>% keep(~  !(. %in% c("eff_v", "MaxE", "RMSE", "model", "type", "loglik")))

# compute model nesting and likelihood ratio tests
model_nestings <- rbind(computeRatioTests(liks_stoch, param_est_names_stoch, "stochastic"),
                        computeRatioTests(liks_deter, param_est_names_deter, "deterministic"))


# build nesting tree and compute likelihood ratio tests


library(tidygraph)
library(ggraph)

edges <- filter(model_nestings, isnested == 1) %>% 
  rename(from = model_ref,
         to = model_comp) %>% 
  mutate(from = as.numeric(from),
         from = ifelse(type == "stochastic", from + 8, from),
         to = as.numeric(to),
         to = ifelse(type == "stochastic", to + 8, to))

best_liks <-rbind(liks_stoch %>% select(loglik, loglik.se, model, type),
                  liks_deter %>% select(loglik, model, type) %>% mutate(loglik.se = 0)) %>% 
  group_by(type, model) %>% 
  arrange(desc(loglik)) %>% 
  slice(1)

nodes <-  best_liks %>% 
  ungroup %>% 
  mutate(id = as.numeric(model)) %>% 
  select(id, model, loglik, loglik.se, type) %>%
  mutate(id = ifelse(type == "stochastic", id +8, id))

nodes$x <- rep(c(1, 1.75, .75, 1.25, 1, 1.5, 2, 1.75), times = 2)

nodes$nparam <- foreach(r = iter(best_liks, by = "row"), .combine = c) %do% {
  if(r[["type"]] == "deterministic") {
    x <- filter(liks_deter, loglik == r$loglik, model == r$model) %>%  slice(1)
    sum(x[param_est_names_deter]>0)
  } else {
    x <- filter(liks_stoch, loglik == r$loglik, model == r$model) %>% slice(1)
    sum(x[param_est_names_stoch]>0)
  }
}


nodes$y <-nodes$nparam
model_network <- tbl_graph(edges = edges,
                           nodes = nodes,
                           directed = T) 


# plot the plot
pg_labels <- tribble(
  ~type, ~label,
  "deterministic", "A",
  "stochastic", "B"
) %>% 
  mutate(x = min(nodes$x), y = max(nodes$y) + 0.5)

deltax <- 0.1

p.graph <- ggraph(model_network) +
  # geom_node_point(aes(x = x, y = y, size = nparam, color = loglik)) + 
  geom_edge_link(aes(color = LLtest2 , linetype = !LLtest2), width = .8) +
  geom_node_label(aes(label = str_c(model, sprintf("%.2f",loglik), sep = "\n")), size = 3.2) +
  scale_edge_color_manual(values = c("grey", "black"), guide = "none") +
  scale_edge_linetype(guide = "none") +
  scale_y_continuous(breaks = seq(7, 13), limits = c(6.8, 13.5)) +
  scale_x_continuous(limits = c(min(nodes$x) - deltax, max(nodes$x) + deltax)) +
  geom_text(data = pg_labels, aes(x = x, y = y, label = label), size = 7) +
  facet_wrap(~type) +
  theme_graph() +
  labs(x = "" , y = "# of paramters") +
  theme(panel.grid.major.y = element_line(color = "#DED8D8"),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_blank()) 

p.graph
ggsave(p.graph, file = "results/figures/model_nesting_LLratio_test.png", width = 9, height = 5, dpi = 300)



# Compare outputs ---------------------------------------------------------

# DETERMINISTIC MODEL
ode_cases <- read_csv("data/ODE_dates.csv", col_names = F) %>% 
  set_colnames("date") %>% 
  mutate(date = date - 1) %>% 
  bind_cols(read_csv("results/ODE_trajectory_MN.csv", col_names = F)[101,] %>% 
              unlist() %>% 
              tibble(cases = .)) 


doMC::registerDoMC(8)

sim_deter_quantiles <- foreach(model = c("MN", "ME", "MEC"), 
                               .combine = rbind) %do% {
                                 # load trajectories
                                 ode_out <- read_csv(str_c("results/ODE_trajectory_", model ,".csv"), col_names = F)
                                 # remove last row containing the data
                                 ode_out <- ode_out[-nrow(ode_out), ]
                                 
                                 # run simulations from the deterministic models
                                 nsim_deter <- 100
                                 
                                 sim_deter <- foreach(it = iter(ode_out, by = "row"), 
                                                      ic = icount(),
                                                      .combine = cbind) %dopar% {
                                                        matrix(rpois(nrow(ode_cases) * nsim_deter, unlist(it)), ncol = nsim_deter) %>% 
                                                          as.data.frame() %>% 
                                                          set_colnames(str_c("sim", ic, "_", seq(nsim_deter)))
                                                      } %>% 
                                   as_tibble()
                                 
                                 foreach(it = iter(sim_deter, by = "row"), 
                                         .combine = rbind) %dopar% {
                                           x <- unlist(it)
                                           tibble(
                                             q05 = quantile(x, 0.025, na.rm = T),
                                             mean = mean(x, na.rm = T),
                                             q50 = quantile(x, 0.5, na.rm = T),
                                             q95 = quantile(x, 0.975, na.rm = T)
                                           )
                                         } %>% 
                                   bind_cols(ode_cases) %>% 
                                   mutate(model = model)
                               } %>% 
  mutate(type = "deterministic")

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




# get MLE paramter sets for the SIRB and SIRB-Exposure models for final plot 
best_param <- best_liks %>% 
  filter(type == "stochastic", model %in% c("MN", "ME", "MEC")) %>% 
  arrange(desc(loglik)) %>% 
  ungroup %>% 
  left_join(liks_stoch) %>% 
  select(model, one_of(names(coef(sirb_cholera)))) 

# This is to avoid the problem of the rainfall normalization
rain <- read_csv("data/rainfall_data.csv") %>% 
  mutate(date = as.Date(date, format = "%d-%b-%y"),
         time = dateToYears(date),
         rain_std = rain/max(rain))

# re-create pomp object for simulations
sirb_cholera_maxrain <- pomp(sirb_cholera,
                             covar = rain %>%
                               filter(time >= min(sirb_cholera@tcovar) & time <= max(sirb_cholera@tcovar)) %>% 
                               select(time, rain_std) %>%
                               rename(rain = rain_std),
                             tcovar = "time")

# run simulations for each model
sim_stochastic <- foreach(r = iter(best_param, by = "row"), 
                          .combine = rbind) %dopar% {
                            
                            if(r$model == "ME") {
                              simulatePOMP(sirb_cholera, unlist(r[names(coef(sirb_cholera))]), 10000) %>% 
                                mutate(model = r$model)
                            } else {
                              simulatePOMP(sirb_cholera_maxrain, unlist(r[names(coef(sirb_cholera))]), 10000) %>% 
                                mutate(model = r$model)
                            }
                            
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

psim_labels <- tribble(
  ~model, ~type, ~label,
  "MN", "deterministic", "A",
  "MN", "stochastic", "B",
  "ME", "deterministic", "C",
  "ME", "stochastic", "D",
  "MEC", "deterministic", "E",
  "MEC", "stochastic", "F"
) %>% 
  mutate(date = as.Date("2015-06-09"), y = 115,
         model = factor(model, levels = c("MN", "ME", "MEC")))

p.sim <- ggplot(data = rbind(sim_stochastic_quantiles, sim_deter_quantiles),
                aes(x = date))+
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.1, color = simcol, fill = simcol) +
  geom_line(aes(y = q50), color = simcol) +
  geom_line(aes(y = mean), linetype = 2, color = simcol) +
  geom_line(aes(y = cases), color = datacol, lwd = 0.2) +
  geom_point(aes(y = cases), color = datacol, size = 0.8) +
  geom_text(data = psim_labels, aes (y = y, label = label), size = 7) +
  facet_grid(model~type) +
  scale_x_date(date_labels = "%b-%y", expand = c(0,0), limits = as.Date(c("2015-06-01", "2015-09-28"))) +
  scale_y_continuous(expand = c(0,0))+ 
  labs(y = "daily cholera cases", x = "date") +
  theme(panel.grid.major = element_line(color = "lightgray"),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        strip.text = element_blank(),
        axis.title = element_text())

p.sim
ggsave(p.sim, filename = "results/figures/simulations_comparison.png", width = 9, height = 10, dpi = 300)


# Write BICs to file
summary_stoch <- nodes %>% 
  mutate(BIC = -2 * loglik + log(length(sirb_cholera_maxrain@data)) * nparam,
         BF = 1/exp(0.5 * (BIC - min(BIC)))) %>% 
  filter(type == "stochastic") %>% 
  select(model, loglik, loglik.se, nparam, BIC, BF) 

write_csv(summary_stoch, path = "results/stochastic_models_summary.csv")


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


# Sumary statistics models ------------------------------------------------

library(xtable)
model_summary <- summary_stoch %>% 
  inner_join(read_csv(file = "results/deter_models_summary.csv") %>% 
               mutate(loglik.se = 0) %>% 
               select(one_of(colnames(summary_stoch)))
             , by = "model", ., suffix = c(".deter", ".stoch"))

print(xtable(model_summary, 
             align = c("c","l", rep("c", 10)),
             digits = c(1,1, rep(c(2,3,2,2,-1), times = 2))),
      file = "results/model_summary_latex.txt")
