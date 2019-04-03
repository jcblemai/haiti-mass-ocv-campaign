from rpy2 import robjects
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.dates as mdates
import pandas as pd
import numpy as np
import datetime
import yaml
import os

# Convert pandas dataframe
from rpy2.robjects import pandas2ri
pandas2ri.activate()

output_dir = 'output_12-20-2gammaMOD/'
dept_avail = os.listdir(output_dir)


rainfall = pd.read_csv('haiti-data/fromAzman/rainfall.csv', index_col = 0, parse_dates = True)
cases    = pd.read_csv('haiti-data/fromAzman/cases_corrected.csv', index_col=0, parse_dates =True)


compartments = ['S', 'I', 'A', 'RA1', 'RA2', 'RA3', 'RI1', 'RI2', 'RI3', 'W', 'B', 'cases', 'C',
                 "VSd", "VRI1d", "VRI2d", "VRI3d", "VRA1d", "VRA2d", "VRA3d",
                 "VSdd", "VRI1dd", "VRI2dd", "VRI3dd", "VRA1dd", "VRA2dd", "VRA3dd",
                 "VSd_alt", "VRI1d_alt", "VRI2d_alt", "VRI3d_alt", "VRA1d_alt", "VRA2d_alt", "VRA3d_alt",
                 "VSdd_alt", "VRI1dd_alt", "VRI2dd_alt", "VRI3dd_alt", "VRA1dd_alt", "VRA2dd_alt", "VRA3dd_alt"]

stream = open('haiti-data/input_parameters.yaml', 'r')
input_parameters = yaml.load(stream)

dept_name = [list(pop.keys())[0] for pop in input_parameters['population']]

t_start = input_parameters['t_start']
t_for = datetime.date(2029,12,21)
nsim = 10
run_lvl = 4 


class VaccinationScenario():

    def __init__(self, course_year, percent_completely_unvaccinated, percent_onedose, percent_twodoses, not_dep = []):
        pop = {'Artibonite':1727524,
        'Centre':746236,
        'Grande_Anse':468301,
        'Nippes':342525,
        'Nord':1067177,
        'Nord-Est':393967,
        'Nord-Ouest':728807,
        'Ouest':4029705,
        'Sud':774976,
        'Sud-Est':632601}

        ocv_order = ['Centre', 'Artibonite','Ouest','Nord-Ouest','Nord','Sud', 'Nippes','Nord-Est', 'Sud-Est','Grande_Anse']


        self.t_vacc_start = {}
        self.t_vacc_end = {}
        self.p1d_reg = {}
        self.r_v_year = {}

        #20% completely unvaccinated, 10% one-dose only, 70% two doses

        t_init = datetime.date(2018,7,14)
        days_per_departement = int((course_year*365)/len(ocv_order))

        for i, dp in enumerate(ocv_order):
            if dp not in not_dep:
                self.t_vacc_start[dp] = t_init + datetime.timedelta(days=i*days_per_departement)
                self.t_vacc_end[dp] = t_init + datetime.timedelta(days=(i+1)*days_per_departement)
                self.p1d_reg[dp] = percent_onedose/(percent_onedose + percent_twodoses)
                self.r_v_year[dp] = pop[dp]*(100-percent_completely_unvaccinated)/100/days_per_departement * 365.25

            else:
                self.t_vacc_start[dp] = t_init + datetime.timedelta(days=i*days_per_departement)
                self.t_vacc_end[dp] = t_init + datetime.timedelta(days=(i+1)*days_per_departement)
                self.p1d_reg[dp] = 0
                self.r_v_year[dp] =0



S1 = VaccinationScenario(2, 20, 10, 70)
S2 = VaccinationScenario(2, 20, 10, 70, not_dep=['Ouest','Nord-Ouest','Sud', 'Nippes','Nord-Est', 'Sud-Est','Grande_Anse'])
S3 = VaccinationScenario(5, 20, 10, 70)
#S2 = VaccinationScenario(2, 40, 20, 40, not_dep=['Ouest','Nord-Ouest','Sud', 'Nippes','Nord-Est', 'Sud-Est','Grande_Anse'])
S4 = VaccinationScenario(2, 3.33, 1.67, 95)
S5 = VaccinationScenario(2, 3.33, 1.67, 95, not_dep=['Ouest','Nord-Ouest','Sud', 'Nippes','Nord-Est', 'Sud-Est','Grande_Anse'])
S6 = VaccinationScenario(5, 3.33, 1.67, 95)
scenario = S6
scenario_str = 'S6'

n =  10


# Running simulation:
index = pd.DatetimeIndex(start =  t_start,  end = t_for, freq = 'W-SAT')
r_source = robjects.r['source'];

all_sim =  pd.DataFrame(0, index = pd.DatetimeIndex(start =  t_start, 
                        end = t_for, freq = 'W-SAT'), columns = np.arange(10))

vacc_init = pd.read_csv('vacc_init' + scenario_str +  '.csv', index_col = 0, parse_dates = True)
vacc_init.to_csv('covar_mob.csv', index_label='date')
for sim in range(n):
    all_data_vacc = {}
    csv_all_q50 =  pd.DataFrame(0, index = pd.DatetimeIndex(start =  t_start, 
                                      end = t_for, freq = 'W-SAT'), columns = dept_name)

    for i, dp in enumerate(dept_avail):
        dept_data = {}
        robjects.r('departement <- "' + dp + '"')
        robjects.r('output_dir <- "' + output_dir + '"')
        robjects.r('run_level <- ' + str(run_lvl))
        robjects.r('nsim <- ' + str(nsim))
    
        robjects.r('t_vacc_start <- "' + str(scenario.t_vacc_start[dp]) + '"')
        robjects.r('t_vacc_end  <- "' + str(scenario.t_vacc_end[dp]) + '"')
        robjects.r('p1d_reg <- ' + str(scenario.p1d_reg[dp]))
        robjects.r('r_v_year <- ' + str(scenario.r_v_year[dp]))
        robjects.r('cases_ext <- 1')
        if dp == 'Artibonite':
            robjects.r('calib_corr <- ' + str((cases[t_start:].mean().sum() - cases[t_start:][dp].mean())/10))
        elif (dp == 'Nord-Est'):
            robjects.r('calib_corr <- ' + str((cases[t_start:].mean().sum() - cases[t_start:][dp].mean())*8))
        elif ((dp == 'Sud-Est') or (dp == 'Nippes')):
            robjects.r('calib_corr <- ' + str((cases[t_start:].mean().sum() - cases[t_start:][dp].mean())*5))
        elif dp == 'Ouest':
            robjects.r('calib_corr <- ' + str((cases[t_start:].mean().sum() - cases[t_start:][dp].mean())/15))
        else:
            robjects.r('calib_corr <- ' + str(cases[t_start:].mean().sum() - cases[t_start:][dp].mean()))

        print(sim, dp)

        r_source('scripts/forecast_haitiOCV_mob.R')
  
        for comp in compartments:
            temp = pandas2ri.ri2py(robjects.r[comp])
            temp.index = index
            temp.drop('date',axis=1, inplace = True)
            dept_data[comp] = temp
        all_data_vacc[dp] = dept_data
        csv_all_q50[dp] = all_data_vacc[dp]['C']['q50']
    all_sim[sim] = csv_all_q50.sum(axis=1)
        
all_sim.to_csv('all_sim' + scenario_str +  '.csv', index_label='date')        
