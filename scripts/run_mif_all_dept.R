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
print(args[1])
run_level <- as.integer(args[1])

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
load(paste0(output_dir, "sirb_cholera_pomped_all.rda"))

# run on a single multi-core machine or on a cluster (using SLURM)
ncpus <- detectCores()
jobname <- "HaitiOCVAll"
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

# values of the random walks standard deviations [From somewhere, don't touch]
rw.sd_rp <- 0.02  # for the regular (process and measurement model) paramters
rw.sd_ivp <- 0.2  # for the initial value paramters
rw.sd_param <- set_names(c(rw.sd_rp, rw.sd_ivp), c("regular", "ivp"))

# Level of detail on which to run the computations [Allow to chose easly set of params]
cholera_Np <-           c(3e3,    3e3,    10e3,    4e3)
cholera_Nmif <-         c(300,     300,   100,    150)      # Entre 200 et 300  
cholera_Ninit_param <-  c(n_runs,  n_runs, n_runs, n_runs)   # How many rounds a cpu does
cholera_NpLL <-         c(1e4,    1e4,    1e4,    1e4)      # Au moins 10 000 pour un truc ok
cholera_Nreps_global <- c(5,      10,      10,     10)


# Run the computations -----------------------------------------------
print(sprintf(">>>> Running on all departements with run level %d", run_level))
print(sprintf(">>>> Np : %d | Nmif: %d | Ninit: %d | NpLL: %d | Nrep: %d ", 
              cholera_Np[run_level], cholera_Nmif[run_level], cholera_Ninit_param[run_level], 
              cholera_NpLL[run_level], cholera_Nreps_global[run_level]))


# names of results files
mifruns.filename = str_c(output_dir, str_c(projname, run_level, sep = "-"), "-mif_runs.rda", sep = "")

# create random vectors of initial paramters given the bounds
best_param <- read_csv("starting.csv") 
to_generate = cholera_Ninit_param[run_level] - nrow(best_param)

allready_there = nrow(best_param)
  
if (to_generate > 0) {
    for (i in 0:(to_generate-1)) {
      new_par = as.data.frame(t(rtruncnorm(ncol(best_param), 
                                           a=0, 
                                           b=Inf,
                                           mean = unlist(slice(best_param,  i%%allready_there +1 )), 
                                           sd = unlist(slice(best_param,  i%%allready_there  +1)/10))))  
      names(new_par) <- names(best_param)
      best_param <- bind_rows(best_param, new_par)
    }
} else if (to_generate < 0) {
    best_param <- head(best_param, to_generate)
}

  
init_params <- as.data.frame(best_param)

job_rw.sd <- eval(
  parse(
    text = str_c("rw.sd(",
                 "betaBArtibonite         = ",   rw.sd_param["regular"],
                 ",betaBSud_Est           = ",   rw.sd_param["regular"],
                 ",betaBNippes            = ",   rw.sd_param["regular"],
                 ",betaBNord_Est          = ",   rw.sd_param["regular"],
                 ",betaBOuest             = ",   rw.sd_param["regular"],
                 ",betaBCentre            = ",   rw.sd_param["regular"],
                 ",betaBNord              = ",   rw.sd_param["regular"],
                 ",betaBSud               = ",   rw.sd_param["regular"],
                 ",betaBNord_Ouest        = ",   rw.sd_param["regular"],
                 ",betaBGrande_Anse       = ",   rw.sd_param["regular"],
                 ",foi_addArtibonite      = ",   rw.sd_param["regular"],
                 ",foi_addSud_Est         = ",   rw.sd_param["regular"],
                 ",foi_addNippes          = ",   rw.sd_param["regular"],
                 ",foi_addNord_Est        = ",   rw.sd_param["regular"],
                 ",foi_addOuest           = ",   rw.sd_param["regular"],
                 ",foi_addCentre          = ",   rw.sd_param["regular"],
                 ",foi_addNord            = ",   rw.sd_param["regular"],
                 ",foi_addSud             = ",   rw.sd_param["regular"],    
                 ",foi_addNord_Ouest      = ",   rw.sd_param["regular"],
                 ",foi_addGrande_Anse     = ",   rw.sd_param["regular"],
                 ", mu_B   = ", 1/10*  rw.sd_param["regular"],
                 ", XthetaA= ", 1/10*  rw.sd_param["regular"],
                 ", thetaI = ", 1/10*  rw.sd_param["regular"],
                 ", lambdaR = ",1/10*  rw.sd_param["regular"],
                 ", r      = ",1/10*   rw.sd_param["regular"],
                 ", std_W  = ",1/10*   rw.sd_param["regular"],
                 ", epsilon= ",1/10*   rw.sd_param["regular"],
                 ", k = "     ,1/10*   rw.sd_param["regular"],
                 ", cas_def = ifelse(time<2018., 0, ",  rw.sd_param["regular"], ")",
                 ")")
  )
)

# get the names of the paramters that are fixed
param_fixed_names <- setdiff(names(coef(sirb_cholera)), colnames(init_params))

# bind with the fixed valued paramters
init_params <- cbind(init_params, 
                     matrix(rep(coef(sirb_cholera)[param_fixed_names],
                                each = cholera_Ninit_param[run_level]),
                            nrow = cholera_Ninit_param[run_level]) %>% 
                       set_colnames(param_fixed_names))


init_param_tocompute <- init_params



# Run MIF
tic("MIF")
# file to store all explorations of the likelihood surface
all_loglik.filename <- sprintf("%sHaiti_OCVAll-param_logliks-10-l%i.csv", output_dir, run_level)
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
           cooling.fraction.50 = 0.5,   # 0.4-> Stabilize after 200 -> stop at 300.
           transform = TRUE,
           rw.sd = job_rw.sd,
           verbose = T
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

closeAllConnections()

# Send a telegram message

if (hostname != 'echopc27') {
  system(sprintf('bash ~/science_bot.sh "Successful completion on %s :) \n (Departement: All | Run level: %d)"', hostname, run_level))
}



