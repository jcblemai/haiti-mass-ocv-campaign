# Rainfall effects on the 2015 Juba cholera epidemic


## TODO
  - *Demographic stochasticity*
  - Beta Drop to match the drop with two beta
  - Use artibonite (Xcept rainfall) as starting point for other simulation
  - Rain may does not account for mobility -> seems ok
  
  
  1. Implement mobility 30 minscenearios
  2. Create good param Arti. Relaunch calib on other from that 1h. vary beta and lambda
  3. Create code to generate 
  4. *VE proportional to population U5* 
  5. Wake up early and generate results.



### TODO vaccination
  - epsilon might be 30%
  +  - Why isn't the first value of the dataframe the same between q05 and q95

  

# As of the last big calibration we have (12-20) is:
  - good in Artibonite, Sud, Grande Anse, Centre, Sud-Est
  - bad in Nord-Est, Nord, Mord Ouest, OUest, Centre, Nippes.

Drop is needed in Nippes, Grande Anse, Nord, Nord-Est, Ouest, Sud, ...


## Data

1. `cases_data.csv`: reported cholera cases for 2014-2015 in Juba
2. `rain_data.csv`: remote sensing estimate of daily rainfall (TRMM data)
3. `input_parameters.yaml`: fixed parameters for the SIRB model in [YAML](http://yaml.org/) format

## Code

1. `pomp_cholera_juba.R`: build the [POMP](https://kingaa.github.io/pomp/) object containing data and code for simulation
2. `run_mif_cholera.R`: fit the models using [multiple iterated filtering](http://www.pnas.org/content/112/3/719)
3. `analysis_mif_runs.R`: code to produce all the figures of the paper, maybe the most usefull is the one to simulate the stochastic model in the "Compare outputs" section


## Package to install

  install.packages(c("tictoc", "pomp", "tidyverse", "magrittr", "ggthemes", "GGally", "foreach", "itertools", "lubridate", "dplyr", "purrr", "readr", "stringr", "tibble", "doMC","doSNOW"  , dependencies=TRUE)

## Show error in R

  options(show.error.locations=TRUE)
  options(error=recover)
