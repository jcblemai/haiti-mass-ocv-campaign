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
rm(list = ls())
Sys.setlocale("LC_ALL","C")

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

# simple function to easily set values of paramters to retain
setRandomWalkSD <- function(job_param_specs, param, param.type = "regular", rw.sd = rw.sd_param) {
  ifelse(job_param_specs[param], rw.sd_param[param.type], 0)
}

# Load pomp object ---------------------------------------------------------------
cholera_pomp_file <- "data/sirb_cholera_pomped.rda"

if(!file.exists(cholera_pomp_file)) {
  # if not done yet create the pomp object
  source("scripts/pomp_cholera_juba.R")
} else {
  load(cholera_pomp_file)
}


# Define model structures -------------------------------------------------
# abbreviations of process names
h2h <- "human-to-human" 
bacteria <- "environmental bacteria reservoir"
rain_expo <- "rainfall effect on exposure"
rain_cont <- "rainfall effect on contamination"

# Specify model types and paramters to estimate
model_specs <- tribble(
  ~model, ~processes, 
  "SIR_B", bacteria,
  "SIR_HB", str_c(h2h, bacteria, sep = " + "),
  "SIR_B_E", str_c(bacteria, rain_expo, sep = " + "),
  "SIR_B_R", str_c(bacteria, rain_cont, sep = " + "),
  "SIR_B_ER", str_c(bacteria, rain_expo, rain_cont, sep = " + "),
  "SIR_HB_E", str_c(h2h, bacteria, rain_expo, sep = " + "),
  "SIR_HB_R", str_c(h2h, bacteria, rain_cont, sep = " + "),
  "SIR_HB_ER", str_c(h2h, bacteria, rain_expo, rain_cont, sep = " + ")
)

# Speficy the paramter names corresponding to the model paramters
param_specs <- tribble(
  ~process, ~param_name,
  h2h, "beta_I",
  bacteria, "beta_B",
  bacteria, "mu_B",
  bacteria, "theta",
  rain_expo, "lambda_E",
  rain_expo, "alpha_E",
  rain_cont, "lambda_R",
  rain_cont, "alpha_R"
)

# Specify which parameters to estimate for each model
for(r in seq(nrow(param_specs))) {
  model_specs[[param_specs$param_name[r]]] <- map_lgl(model_specs$processes, 
                                                      ~param_specs$process[r] %in% str_split(., " \\+ ")[[1]])
}

# Parallel setup ----------------------------------------------------------

# run on a single multi-core machine or on a cluster (using SLURM)
oneMACHINE  = T

if (!oneMACHINE) {
  jobname <- Sys.getenv("SLURM_JOB_NAME")
  job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")
  runid <- Sys.getenv("SLURM_JOB_ID")
  ntasks <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
  ncpus <- as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  rscript.options = commandArgs()
  script = sub(".R","",sub(".* = ", "",rscript.options[grep("--file",rscript.options)]))
  
  # seed for simulations to assure reproducibility
  master.seed <- as.integer(runid) 
  
  # the array id is used to specify the model type to run
  array_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  array_id_vec <- array_id
  # spawn workers
  registerDoMC(ncpus)
  
  n_runs_pertask <- ncpus * 3
  
} else {
  ncpus <- detectCores()
  jobname <- "choleraJuba"
  script <- "run_mif_cholera"
  projname <- jobname
  job_id <- 1
  # seed for simulations to assure reproducibility
  master.seed <- as.integer(runif(1) * 10000) 
  
  # number of runs to do for each task (n * number of cores)
  n_runs_pertask <- ncpus * 1
  
  # simple setup to run only one initial condition per model per node
  ntasks <-  n_runs_pertask * nrow(model_specs) 
  
  # array of ids to run
  array_id_vec <- seq(1,1) * 100 + 1
  # spawn workers
  registerDoMC(ncpus)
  
}

# Set paramter bounds -----------------------------------------------------

# lower bound for positive parameter values
min_param_val <- 1e-5 
# define the bounds for the paramters to estimate
parameter_bounds <- tribble(
  ~param, ~lower, ~upper,
  # Regular paramters
  "sigma", min_param_val, 1 - min_param_val,
  "beta_B", min_param_val, 50,
  "beta_I", min_param_val, 1e3,
  "mu_B", min_param_val, 1e2,
  "theta", min_param_val, 1e2,
  "lambda_E", min_param_val, 5,
  "lambda_R", min_param_val, 5,
  "alpha_E", min_param_val, 40,
  "alpha_R", min_param_val, 40,
  "rho", 0.02, 20,
  # Process noise
  "std_W", min_param_val, 1e-1,
  # Measurement model
  "epsilon", min_param_val, 2,
  #"k", min_param_val, 10,
  "R_0", min_param_val, 0.2
)

# convert to matrix for ease
parameter_bounds <- set_rownames(as.matrix(parameter_bounds[, -1]), parameter_bounds[["param"]])

# Setup MIF paramters -----------------------------------------------------

# values of the random walks standard deviations
rw.sd_rp <- 0.02  # for the regular (process and measurement model) paramters
rw.sd_ivp <- 0.2  # for the initial value paramters
rw.sd_param <- set_names(c(rw.sd_rp, rw.sd_ivp), c("regular", "ivp"))

# Level of detail on which to run the computations
run_level <- 1
cholera_Np <- c(1000, 3e3, 1e4)
cholera_Nmif <- c(1, 300, 400)
cholera_Ninit_param <- c(n_runs_pertask, n_runs_pertask, 10)
cholera_NpLL <- c(2000, 1e4, 5e4)
cholera_Nreps_global <- c(1, 20, 100)

# Run the computations -----------------------------------------------
# if on one machine run all models in sequence if only one task per model 
# seq(1, nrow(model_specs))
for(array_id in array_id_vec) {
  
  # get model for current job array ID: 
  model_id <- floor(array_id/100)
  task_id <- array_id - model_id * 100 
  
  # number of taskes per model
  ntasks_permodel <- floor(ntasks/nrow(model_specs))
  
  # select model for this job in array
  job_model_specs <- model_specs[model_id, ]
  job_param_specs <- unlist(job_model_specs[, -(1:2)])
  
  # names of results files
  mifruns.filename = str_c("results/", str_c(projname, job_model_specs$model, job_id, array_id, str_c(run_level, sep = "-")), "-mif_runs.rda", sep = "")
  
  # set to 0 the paramters that should not be fit
  parameter_bounds[names(job_param_specs)[!job_param_specs], ] <- 0
  
  # create random vectors of initial paramters given the bounds
  init_params <- sobolDesign(lower = parameter_bounds[, "lower"],
                             upper = parameter_bounds[, "upper"], 
                             nseq = ntasks_permodel * cholera_Ninit_param[run_level])
  
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
                   ", beta_B  = ",  setRandomWalkSD(job_param_specs, "beta_B"),
                   ", beta_I  = ",  setRandomWalkSD(job_param_specs, "beta_I"),
                   ", mu_B  = ",  setRandomWalkSD(job_param_specs, "mu_B"),
                   ", theta  = ",  setRandomWalkSD(job_param_specs, "theta"),
                   ", rho  = ",  rw.sd_param["regular"],
                   ", lambda_E  = ",  setRandomWalkSD(job_param_specs, "lambda_E"),
                   ", lambda_R  = ",  setRandomWalkSD(job_param_specs, "lambda_R"),
                   ", alpha_E  = ",  setRandomWalkSD(job_param_specs, "alpha_E"),
                   ", alpha_R  = ",  setRandomWalkSD(job_param_specs, "alpha_R"),
                   ", std_W  = ",  rw.sd_param["regular"],
                   ", epsilon  = ",  rw.sd_param["regular"],
                   ", R_0  = ivp(",  rw.sd_param["ivp"], ")",
                   ")")
    )
  )
  
  if(!oneMACHINE){
    # cut up the initial parameter matrices for computations
    init_params_ichunk <- ichunk(iter(init_params, "row"), cholera_Ninit_param[run_level])
    
    # get the task to run in this instance of the array
    init_param_tocompute <- foreach(chnk = init_params_ichunk) %do% {
      bind_rows(chnk)
    }
    
    init_param_tocompute <- init_param_tocompute[[task_id]]
  } else {
    init_param_tocompute <- init_params
  }
  
  # Run MIF
  # file to store all explorations of the likelihood surface
  all_loglik.filename <- sprintf("results/choleraJuba_param_logliks-10-l%i.csv", run_level)
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
    
    # Compute likelihood
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
    } %>%
      mutate(model = job_model_specs$model)
    
    # write the results to a global file with all preceeding runs for each run level
    write_csv(lik_mf, all_loglik.filename, append = file.exists(all_loglik.filename))
  }, 
  seed = master.seed, kind="L'Ecuyer")
  
}

closeAllConnections()

