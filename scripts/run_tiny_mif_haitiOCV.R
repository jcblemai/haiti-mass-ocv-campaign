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

# TODO manually change:
coef(sirb_cholera)["sigma"] <- 6.651372e-01
coef(sirb_cholera)["mu_B"] <-  8.805273e+02
coef(sirb_cholera)["XthetaA"] <- 3.717466e-01
coef(sirb_cholera)["thetaI"] <- 6.061175e-04
coef(sirb_cholera)["rhoA"] <- 5.051574e+01
coef(sirb_cholera)["XrhoI"] <- 8.996542e-01
coef(sirb_cholera)["std_W"] <- 3.160177e-03
coef(sirb_cholera)["epsilon"] <- 2.706321e-01
coef(sirb_cholera)["k"] <- 8.302246e+01


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
  "betaB", 6.244070e-01, 6.244070e-01,
  "lambdaR", 1.470360e-01, 1.470360e-01,
  "r", 3.913512e-01, 3.913512e-01,
  "foi_add", 7.257634e-07, 7.257634e-07
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
# level 4
cholera_Np <-           c(1e2,    1e3,    2e3,    3e3)
cholera_Nmif <-         c(5,      100,    200,    400)      # Entre 200 et 300  
cholera_Ninit_param <-  c(n_runs, n_runs, n_runs*2, n_runs*4)   # How many rounds a cpu does
cholera_NpLL <-         c(1e2,    1e4,    1e4,    1e4)      # Au moins 10 000 pour un truc ok
cholera_Nreps_global <- c(1,      5,      10,     15)


# Run the computations -----------------------------------------------
# if on one machine run all models in sequence if only one task per model 
# seq(1, nrow(model_specs))
for(array_id in array_id_vec) {
  print(sprintf(">>>> Running on departement %s with run level %d", departement, run_level))
  print(sprintf(">>>> Np : %d | Nmif: %d | Ninit: %d | NpLL: %d | Nrep: %d ", 
                cholera_Np[run_level], cholera_Nmif[run_level], cholera_Ninit_param[run_level], 
                cholera_NpLL[run_level], cholera_Nreps_global[run_level]))
  
  # get model for current job array ID: 
  
  # select model for this job in array
  
  # names of results files
  mifruns.filename = str_c(output_dir, departement, "/", str_c(projname, run_level, departement, sep = "-"), "-mif_runs.rda", sep = "")
  
  # create random vectors of initial paramters given the bounds
  init_params <- sobolDesign(lower = parameter_bounds[, "lower"],
                             upper = parameter_bounds[, "upper"], 
                             nseq = cholera_Ninit_param[run_level])
  
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
      text = str_c("rw.sd(betaB     = ",  rw.sd_param["regular"],
                   ", lambdaR = ",  rw.sd_param["regular"],
                   ", r       = ",  rw.sd_param["regular"],
                   ", foi_add = ",  rw.sd_param["regular"],
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
  
}

closeAllConnections()

# Send a telegram message

if (hostname != 'echopc27') {
  system(sprintf('bash ~/science_bot.sh "Successful completion on %s :) \n (Departement: %s | Run level: %d)"', hostname, departement, run_level))
}



