# Title: Run mif on the Juba cholera model
# Description: Parallel runs of mif2 on different versions of  the model refining first round of results
# Date: Tue Jul 10 18:02:55 2018
# Author: javier.perezsaez@epfl.ch

# Preamble ---------------------------------------------------------------
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(magrittr)
library(tibble)
library(pomp)
library(foreach)
library(itertools)
library(doMC)
library(lubridate)
library(tictoc)
library(truncnorm)
rm(list = ls())
Sys.setlocale("LC_ALL","C")
hostname <- system('hostname', intern = T) 
output_dir <- "output/"

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  # default departement
  args[1] = "Nippes"
  args[2] = 3
} else if (length(args)==1) {
  args[2] = 1
}

# Choose to restart from a previous file (named checkpoint.csv)
restart <- T



departement <- args[1]
run_level <- as.integer(args[2])


# Helper functions 
dateToYears <- function(date, origin = as.Date("2014-01-01"), yr_offset = 2014) {
  julian(date, origin = origin)/365.25 + yr_offset
}

yearsToDate <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
  as.Date((year_frac - yr_offset) * 365.25, origin = origin)
}

yearsToDateTime <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
  as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
}

# Load pomp object ---------------------------------------------------------------
load(paste0(output_dir, departement, "/sirb_cholera_pomped_", departement, ".rda"))

# Parallel setup ----------------------------------------------------------

# run on a single multi-core machine or on a cluster (using SLURM)

ncpus <- detectCores()
jobname <- "HaitiOCV"
script <- "run_mif_cholera"
projname <- jobname
job_id <- 1
# seed for simulations to assure reproducibility
master.seed <- as.integer(runif(1) * 10000) 

# number of runs to do for each task (n * number of cores)
n_runs <- ncpus * 1

# array of ids to run
array_id_vec <- seq(1,1) * 100 + 1
# spawn workers
registerDoMC(ncpus)


# Set parameter bounds -----------------------------------------------------

# lower bound for positive parameter values
min_param_val <- 1e-5 
# define the bounds for the paramters to estimate, juste to give initial parameters.
parameter_bounds <- tribble(
  ~param, ~lower, ~upper,
  # Regular paramters
  "sigma", 0.3, 1 - min_param_val,
  "betaB", min_param_val, 3,
  "mu_B", min_param_val, 1e2,
  "XthetaA", min_param_val, .5,
  "thetaI", min_param_val, 2,
#  "lambda", min_param_val, 5,
   "lambdaR", min_param_val, 5,
  # "gammaA", 73, 365,
  # "gammaI", 73, 365,
  "r", min_param_val, 2,
  #"rhoA", 0.02, 10,
  #"XrhoI", min_param_val, 1,
  # Process noise
  "std_W", min_param_val, 1e-1,
  # Measurement model
  "epsilon", min_param_val, 1,
  "foi_add", min_param_val, 0.005,
  "k", -3, 4 #,   # hard to get negbin like this, sobol in log scale -5 et 4
#  "Rtot_0", min_param_val, 0.1
)


# convert to matrix for ease
parameter_bounds <- set_rownames(as.matrix(parameter_bounds[, -1]), parameter_bounds[["param"]])

# Setup MIF paramters -----------------------------------------------------

# values of the random walks standard deviations [From somewhere, don't touch]
rw.sd_rp <- 0.02  # for the regular (process and measurement model) paramters
rw.sd_ivp <- 0.2  # for the initial value paramters
rw.sd_param <- set_names(c(rw.sd_rp, rw.sd_ivp), c("regular", "ivp"))

# Level of detail on which to run the computations [Allow to chose easly set of params]
# level 1 is short
# level 2 is 12h
# level 4 is 3.2 days
cholera_Np <-           c(1e2,    3e3,    3e3,    4e3)
cholera_Nmif <-         c(5,      300,    300,    400)      # Entre 200 et 300  
cholera_Ninit_param <-  c(n_runs, n_runs, n_runs*2, n_runs*4)   # How many rounds a cpu does
cholera_NpLL <-         c(1e2,    1e4,    1e4,    1e4)      # Au moins 10 000 pour un truc ok
cholera_Nreps_global <- c(1,      5,      10,     15)


# Run the computations -----------------------------------------------
# if on one machine run all models in sequence if only one task per model 
# seq(1, nrow(model_specs))
#for(array_id in array_id_vec) {
  print(sprintf(">>>> Running on departement %s with run level %d", departement, run_level))
  print(sprintf(">>>> Np : %d | Nmif: %d | Ninit: %d | NpLL: %d | Nrep: %d ", 
                cholera_Np[run_level], cholera_Nmif[run_level], cholera_Ninit_param[run_level], 
                cholera_NpLL[run_level], cholera_Nreps_global[run_level]))
  
  # get model for current job array ID: 
  
  # select model for this job in array

  # names of results files
  mifruns.filename = str_c(output_dir, departement, "/", str_c(projname, run_level, departement, sep = "-"), "-mif_runs.rda", sep = "")
  
  # create random vectors of initial paramters given the bounds
  if (restart) {
    liks_stoch <- read_csv("checkpoint.csv") 
    best_like <- liks_stoch %>% 
           arrange(desc(loglik)) %>% 
           slice(1)
    best_like <- best_like$loglik

    # get MLE paramter sets 
    best_param <- liks_stoch %>% filter(loglik > best_like - 2) %>% arrange(desc(loglik)) 
    
    to_generate = cholera_Ninit_param[run_level] - nrow(best_param)
    
    print(to_generate)
    allready_there = nrow(best_param)
    
    if (to_generate > 0) {
      for (i in 0:(to_generate-1)) {
        new_par = as.data.frame(t(rtruncnorm(ncol(best_param), 
                                    a=0, 
                                    b=Inf,
                                    mean = unlist(slice(best_param,  i%%allready_there +1 )), 
                                    sd =   unlist(slice(best_param,  i%%allready_there  +1)/1))))  
        names(new_par) <- names(best_param)
        best_param <- bind_rows(best_param, new_par)
      }
    } else if (to_generate < 0) {
      best_param <- head(best_param, to_generate)
    }
    
    
    best_param <- best_param[rownames(parameter_bounds)]
    
    init_params <- as.data.frame(best_param)
    
  } else {
    init_params <- sobolDesign(lower = parameter_bounds[, "lower"],
                               upper = parameter_bounds[, "upper"], 
                               nseq = cholera_Ninit_param[run_level])
    
    # Allow large variation of k to chose neg in and poisson
    init_params %<>% mutate(k=10^k) 
    
  }

  
  # get the names of the paramters that are fixed
  param_fixed_names <- setdiff(names(coef(sirb_cholera)), colnames(init_params))
  
  # bind with the fixed valued paramters
  init_params <- cbind(init_params, 
                       matrix(rep(coef(sirb_cholera)[param_fixed_names],
                                  each = cholera_Ninit_param[run_level]),
                              nrow = cholera_Ninit_param[run_level]) %>% 
                         set_colnames(param_fixed_names))
  
  
  # set random walk parameters for MIF
  # Define the variance of the perturbation kernel for the paramters
  job_rw.sd <- eval(
    parse(
      text = str_c("rw.sd(sigma  = ",  rw.sd_param["regular"],
                   ", betaB  = ",  rw.sd_param["regular"],
                   ", mu_B   = ",  rw.sd_param["regular"],
                   ", XthetaA= ",  rw.sd_param["regular"],
                   ", thetaI = ",  rw.sd_param["regular"],
                   #", lambda = ",  rw.sd_param["regular"],
                   ", lambdaR = ",  rw.sd_param["regular"],
                   ", r      = ",  rw.sd_param["regular"],
                   #", XrhoI  = ",  rw.sd_param["regular"],
                   #", rhoA   = ",  rw.sd_param["regular"],
                   ", std_W  = ",  rw.sd_param["regular"],
                   #", gammaA  = ", rw.sd_param["regular"],
                   #", gammaI  = ", rw.sd_param["regular"],
                   #", Rtot_0  = ivp(",  rw.sd_param["ivp"],")",
                   ", epsilon= ",  rw.sd_param["regular"],
                   ", foi_add= ",  rw.sd_param["regular"],
                   ", k = "     ,  rw.sd_param["regular"],        # to get binomial, comment for poisson.
                   ")")
    )
  )
  
  init_param_tocompute <- init_params
  
  
  # Run MIF
  tic("MIF")
  # file to store all explorations of the likelihood surface
  all_loglik.filename <- sprintf("%s%s/Haiti_OCV-%s-param_logliks-10-l%i.csv", output_dir, departement, departement, run_level)
  # run computations (stew ensures not to duplicate calculations and sets RNG)
  
  stew(mifruns.filename, {
    w1 <- getDoParWorkers()
    run.seed <- master.seed
    # Run iterated filtering
    t1 <- system.time({
      mf <- foreach(itpar = iter(init_param_tocompute, "row"),
                    .packages = c("pomp"),
                    .combine = c,
                    .errorhandling = "remove",
                    .options.multicore = list(set.seed = TRUE)
      ) %dopar% {
        mif2(sirb_cholera,
             start = unlist(itpar),
             Np = cholera_Np[run_level],
             Nmif = cholera_Nmif[run_level],
             cooling.type = "geometric",
             cooling.fraction.50 = 0.4,
             transform = TRUE,
             rw.sd = job_rw.sd,
             verbose = F
        )
      }
    })
    tt <- toc()
    if (hostname != 'echopc27') {
      system(sprintf('bash ~/science_bot.sh "Successfully done MIF on %s, elapsed time %.2f min :)"', hostname, (tt$toc -tt$tic)/60))
    }
    

    # Compute more precise likelihood (because mif uses a fast like)
    tic("LL")
    lik_mf <- foreach(mifit = mf,
                      i = icount(length(mf)),
                      .packages = c("pomp", "tibble", "dplyr"),
                      .combine = rbind,
                      .errorhandling = "remove",
                      .options.multicore = list(set.seed=TRUE)
    ) %dopar% {
      set.seed(run.seed + i)
      
      ll_mean <- logmeanexp(
        replicate(cholera_Nreps_global[run_level],
                  logLik(
                    pfilter(sirb_cholera,
                            params = pomp::coef(mifit),
                            Np = cholera_NpLL[run_level])
                  )
        ), se = TRUE)
      
      tibble(loglik = ll_mean[1], loglik.se = ll_mean[2]) %>%
        cbind(data.frame(t(coef(mifit))))
    }
    tt <- toc()
    if (hostname != 'echopc27') {
      system(sprintf('bash ~/science_bot.sh "Successfully done LL on %s, elapsed time %.2f min :)"', hostname, (tt$toc -tt$tic)/60))
    }
    
    
    # write the results to a global file with all preceeding runs for each run level
    write_csv(lik_mf, all_loglik.filename, append = file.exists(all_loglik.filename))
  }, 
  seed = master.seed, kind="L'Ecuyer")
  
#}

closeAllConnections()

# Send a telegram message

if (hostname != 'echopc27') {
    system(sprintf('bash ~/science_bot.sh "Successful completion on %s :) \n (Departement: %s | Run level: %d)"', hostname, departement, run_level))
}



