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

output_dir <- 'output/'

load(paste0(output_dir, "/sirb_haitiOCV_pomped.rda"))

sirb_haitiOCV %>%
  simulate(nsim=20,format="data.frame")  -> sims

sims %>%
  ggplot(aes(x=time, y=CasesAll, group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)
