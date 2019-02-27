# Rainfall effects on the 2015 Juba cholera epidemic


## TODO
CODE CALIB
  - Demographic stochasticity ?
  - Beta Drop to match the drop with two beta ?
  - Calibrate Ouest later (using run_tiny_mif) ?


NOTEBOOK AND CODE RUN
  - Mobility a better way:
    1. Manage to recreate initial conditions: launch sim by simulation b.c pomp store state
    2 Do it in the pomp object
    3. Iterative (easier but wronger) via csv
      -> I keep nobility like this because it's too long otherwise
  - *VE proportional to population U5* 
  - Find what's need to be calibrated
  - Why isn't the first value of the dataframe the same between q05 and q95
  - Why all simulation begins at 0 ?

## DONE since feb 20. 2019
  - Use artibonite (Xcept rainfall) as starting point for other simulation
  - add new data
  - Cancel some parameters
  - being able to run stuff in parallel




### FIRST RESULT ROUND FOR HAITI JANUARY
    - data used was 'output_12-20-2gammaMOD/'
    - scripts on echopc40 for probablity of extinction
    - foi_add was scaled :)


# As of the last big calibration we have (12-20) is:
  - good in Artibonite, Sud, Grande Anse, Centre, Sud-Est
  - bad in Nord-Est, Nord, Mord Ouest, OUest, Centre, Nippes.
        -> Drop is needed in Nippes, Grande Anse, Nord, Nord-Est, Ouest, Sud, ...


## Data

1. `cases_corrected.csv`: Cases corrected with some NAs on big drops
2. `rain_data.csv`: remote sensing estimate of daily rainfall (TRMM data)
3. `input_parameters.yaml`: fixed parameters for the SIRB model in [YAML](http://yaml.org/) format

## Code

1. `pomp_cholera_juba.R`: build the [POMP](https://kingaa.github.io/pomp/) object containing data and code for simulation
2. `run_mif_cholera.R`: fit the models using [multiple iterated filtering](http://www.pnas.org/content/112/3/719)
3. `analysis_mif_runs.R`: code to produce all the figures of the paper, maybe the most usefull is the one to simulate the stochastic model in the "Compare outputs" section


## Package to install

  install.packages(c("tictoc", "pomp", "tidyverse", "magrittr", "ggthemes", "GGally", "foreach", "itertools", "lubridate", "dplyr", "purrr", "readr", "stringr", "tibble", "doMC","doSNOW","truncnorm")  , dependencies=TRUE)

## Show error in R

  options(show.error.locations=TRUE)
  options(error=recover)
