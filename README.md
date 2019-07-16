# Haiti Mass Vaccination campaign

## TODO calibration

1. Change FOI
2. Past vaccination campaings ?ß
2. Lanch big calibration


## TODO experiment

1. Tipping point for vaccination
2. Difference rainfall scenerio

### Code: Shell script

1. `run.sh`: run the calibration for a departement and a run level.
1. `generate.sh`: run the calibration or just some generation (like pomp object, parameter file) for all departements at a run level.


### Code: R scripts

1. `pomp_cholera_juba.R`: build the [POMP](https://kingaa.github.io/pomp/) object containing data and code for simulation
2. `run_mif_cholera.R`: fit the models using [multiple iterated filtering](http://www.pnas.org/content/112/3/719)
3. `analysis_haitiOCV.R`: plot some information about calibration.
5. `forecast_haitiOCV.R`: code to project the model. Take files `haiti-data/proj/rainfall.csv` for rainfall and file `covar_mob.csv` for the mobility covariate.
4. `trials.R`: A little scratchpad ! (now to plot also diagnosis)


### Code: Python notebook

1. `forecast.ipynb`: Code to forecast and project vaccination scenarios
2. `rainfall.ipynb`: Create a rainfall projection file, and run some checks if we add new data.
3. `data_analysis.ipynb`: Data analysis on the initial data.


### Code: C functions

1. `sirb_model_vacc.c`: The pomp model containing the transitions, able to support two vaccination campaign.
2. `v_eff.c`: Function for the two doses vaccine efficacity

### Data

In the fromAzman folder:
1. `cases_corrected.csv`: Cases corrected with some NAs on big drops
2. `rainfall.csv`: remote sensing estimate of daily rainfall (TRMM data)
In the proj folder:
1. `rainfall.csv`: used for projection.
At root level:
3. `input_parameters.yaml`: fixed parameters for the SIRB model in [YAML](http://yaml.org/) format


## Useful commands

### Package to install

    install.packages(c("tictoc", "pomp", "tidyverse", "magrittr", "ggthemes", "GGally", "foreach", "itertools", "lubridate", "dplyr", "purrr", "readr", "stringr", "tibble", "doMC","doSNOW","truncnorm", "zoo")  , dependencies=TRUE)
 
### Show errors in R
    
    options(show.error.locations=TRUE)
    options(error=recover)
