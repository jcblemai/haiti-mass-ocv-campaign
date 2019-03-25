# Haiti Mass Vaccination campaign


## TODO

### Calibration

  - Demographic stochasticity ?
  - Beta Drop to match the drop with two beta ?
  - Calibrate Ouest later ?


### Forecast and results
  - Why isn't the first value of the dataframe the same between q05 and q95
  - Why all simulation begins at 0 ?
  - Better mobility ? or at least the ability to see simulation by simulation [WIP in trials.R]
  
### Reading
  - Tidy data

### DONE since feb 20. 2019

  - VE proportional to population U5*
  - Find what's need to be calibrated
  - Use artibonite (Xcept rainfall) as starting point for other simulation
  - add new data
  - Cancel some parameters
  - being able to run stuff in parallel
  - Mobility a better way:
    1. Manage to recreate initial conditions: launch sim by simulation b.c pomp store state
    2 Do it in the pomp object
    3. Iterative (easier but wronger) via csv
      -> I keep nobility like this because it's too long otherwise
  - Build some calibration diagnosis
  - Check Artibonite on echopc20




## Results

### Calibration for next round of results (April 8 2019)
    
  - Big calibration in Artibonite
  - Calibrate only foi and beta in all other departement, take param from artibonite
    
We have good results except in *Ouest* (see `15-05-long`)

### Results for the haiti meeting in January:
    
  - data used was 'output_12-20-2gammaMOD/'
  - scripts on echopc40 for probablity of extinction
  - foi_add was scaled :)
  - see `meeting 27-01-2019 output`
  - Allow run of one simulation for diagnosis
  - Plot national scale.


### As of the last big calibration we have (12-20) is:
  - good in Artibonite, Sud, Grande Anse, Centre, Sud-Est
  - bad in Nord-Est, Nord, Mord Ouest, OUest, Centre, Nippes.
        -> Drop is needed in Nippes, Grande Anse, Nord, Nord-Est, Ouest, Sud, ...

## Files

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

    install.packages(c("tictoc", "pomp", "tidyverse", "magrittr", "ggthemes", "GGally", "foreach", "itertools", "lubridate", "dplyr", "purrr", "readr", "stringr", "tibble", "doMC","doSNOW","truncnorm")  , dependencies=TRUE)
 
### Show errors in R
    
    options(show.error.locations=TRUE)
    options(error=recover)
