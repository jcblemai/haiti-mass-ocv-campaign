from rpy2 import robjects
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.dates as mdates
import multiprocessing as mp
import pandas as pd
import numpy as np
import datetime
import yaml
import os
# Convert pandas dataframe
from rpy2.robjects import pandas2ri
pandas2ri.activate()
# Suppress R warnings in python:
import warnings
from rpy2.rinterface import RRuntimeWarning
warnings.filterwarnings("ignore", category=RRuntimeWarning)

# Warning: Needs an output/Simulations folder

nsim = 1000
n_proc = 18
        
compartments = ['S', 'I', 'A', 'RA1', 'RA2', 'RA3', 'RI1', 'RI2', 'RI3', 'W', 'B', 'cases', 'C',
                 "VSd", "VRI1d", "VRI2d", "VRI3d", "VRA1d", "VRA2d", "VRA3d",
                 "VSdd", "VRI1dd", "VRI2dd", "VRI3dd", "VRA1dd", "VRA2dd", "VRA3dd",
                 "VSd_alt", "VRI1d_alt", "VRI2d_alt", "VRI3d_alt", "VRA1d_alt", "VRA2d_alt", "VRA3d_alt",
                 "VSdd_alt", "VRI1dd_alt", "VRI2dd_alt", "VRI3dd_alt", "VRA1dd_alt", "VRA2dd_alt", "VRA3dd_alt"]

departements = ['Artibonite','Centre','Grande_Anse','Nippes','Nord','Nord-Est','Nord-Ouest','Ouest','Sud','Sud-Est']


# Rainfall generation
rainfall = pd.read_csv('haiti-data/fromAzman/rainfall.csv', index_col = 0, parse_dates = True)
stream = open('haiti-data/input_parameters.yaml', 'r')
input_parameters = yaml.load(stream)
dept_name = [list(pop.keys())[0] for pop in input_parameters['population']]
t_start = input_parameters['t_start']
t_for = datetime.date(2029,12,20)
def project_rain(rainfall, tf):
    nd = 14 #days sampled - must be multiple of 7 d
    dti = rainfall.iloc[0].name.date()
    dtf = rainfall.iloc[-1].name.date()
    rain_prj_index = pd.DatetimeIndex(start =  dtf + datetime.timedelta(1), 
                                      end = tf, freq = 'D')
    rain_prj = np.zeros((rain_prj_index.shape[0], 10))
    # Full years of data available
    years = range(dti.year+1, dtf.year-1)
    # each nd days, assign an al precipitation.
    for i, date in enumerate(pd.date_range(dtf + datetime.timedelta(1), tf, freq = str(nd)+'D')):
        dd = date.day
        if (date.month == 2 and dd == 29):
            dd = 28
        pick = datetime.date(np.random.choice(years), date.month, dd)
        #print(pick, i, rainfall.loc[pd.date_range(pick, pick + datetime.timedelta(nd-1))].values.shape, rain_prj[nd * i: nd * (i+1)].shape)
        rain_prj[nd * i: nd * (i+1)] = rainfall.loc[pd.date_range(pick, pick + datetime.timedelta(nd-1))].values
    rain_prj = pd.DataFrame(rain_prj, index = rain_prj_index, columns = dept_name)    
    return rain_prj
rain_prj = project_rain(rainfall, t_for)
rain = pd.concat((rainfall, rain_prj))
rain.to_csv('haiti-data/proj/rainfall.csv', index_label = 'date')

# Starting the run for all departemens


stream = open('haiti-data/input_parameters.yaml', 'r')
input_parameters = yaml.load(stream)

dept_name = [list(pop.keys())[0] for pop in input_parameters['population']]

t_start = input_parameters['t_start']
t_for = datetime.date(2029,12,20)

index = pd.DatetimeIndex(start =  t_start,  end = t_for, freq = 'W-SAT')

class DeptData():
    def __init__(self):
        self.q05 = pd.DataFrame(0, index = pd.DatetimeIndex(start =  t_start, 
                                      end = t_for, freq = 'W-SAT'), columns = compartments)
        self.q50 = pd.DataFrame(0, index = pd.DatetimeIndex(start =  t_start, 
                                      end = t_for, freq = 'W-SAT'), columns = compartments)
        self.q95 = pd.DataFrame(0, index = pd.DatetimeIndex(start =  t_start, 
                                      end = t_for, freq = 'W-SAT'), columns = compartments)

class VaccinationScenario():
    def __init__(self, course_year, percent_completely_unvaccinated, percent_onedose, percent_twodoses, ve, not_dep = []):
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
        self.ve = ve
        self.not_dep = not_dep

        #20% completely unvaccinated, 10% one-dose only, 70% two doses
        t_init = datetime.date(2019,1,12)
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
def run_sim(scenario_str):
    # see second answer of https://stackoverflow.com/questions/25175530/can-rpy2-code-be-run-in-parallel on why it starts different 
    # R instances

    print ("Working on scenario " + scenario_str)
    all_data = {}
    scenario = scenarios[scenario_str]

    for dp in departements:
        all_data[dp] = DeptData()


    r_source = robjects.r['source']
    dept_data = {}
    r_options = robjects.r['options']
    r_options(warn=-1)
    robjects.r('scenario     <- "' + scenario_str + '"')
    robjects.r('nsim         <- '  + str(nsim))
    robjects.r('t_vacc_start <- list()')
    robjects.r('t_vacc_end <- list()')
    robjects.r('p1d_reg <- list()')
    robjects.r('r_v_year <- list()')


    for dp in departements:
        robjects.r('t_vacc_start${} <- "'.format(dp.replace('-','_')) + str(scenario.t_vacc_start[dp]) + '"')
        robjects.r('t_vacc_end${}   <- "'.format(dp.replace('-','_')) + str(scenario.t_vacc_end[dp]) + '"')
        robjects.r('p1d_reg${}     <- '.format(dp.replace('-','_'))  + str(scenario.p1d_reg[dp]))
        robjects.r('r_v_year${}     <- '.format(dp.replace('-','_'))  + str(scenario.r_v_year[dp]))
    robjects.r('cases_ext    <- '  + str(scenario.ve))
    r_source('scripts/forecast_all_dept.R')
    temp = robjects.r['sim_stochastic']
    
    for dp in departements:
        for comp in compartments:
            all_data[dp].q05[comp] = temp[temp['variable'] == comp + dp.replace('-','_')][temp['isdata']=='simulation']['q05'].values
            all_data[dp].q50[comp] = temp[temp['variable'] == comp + dp.replace('-','_')][temp['isdata']=='simulation']['q50'].values
            all_data[dp].q95[comp] = temp[temp['variable'] == comp + dp.replace('-','_')][temp['isdata']=='simulation']['q95'].values


    all_data['Ouest'].q05['cases'][:datetime.date(2017,6,10)] = np.nan
    all_data['Ouest'].q50['cases'][:datetime.date(2017,6,10)] = np.nan
    all_data['Ouest'].q95['cases'][:datetime.date(2017,6,10)] = np.nan

    print ("Finished run of scenario " + scenario_str)
    
    # Save results
    dir_name = 'output/Results/' + scenario_str + '/'
    try:
        os.makedirs(dir_name)    
        print("Directory " , dir_name ,  " Created ")
    except FileExistsError:
        print("Directory " , dir_name ,  " already exists, writing into.")
    save_result(all_data, scenario_str, folder_name = dir_name)
    
    # Do figure
    ti = input_parameters['t_start']
    tf = t_for

    fig, axes = plt.subplots((len(all_data))//2, 2, figsize=(15,15), squeeze = True, dpi = 200)
    fig.patch.set_facecolor('white')

    axes = axes.flat;
    for i, dp in enumerate(departements):

        axt =  axes[i].twinx()
        axes[i].plot(cases[dp][t_start:][ti:tf], marker='.', linestyle='-',color='black', linewidth=0, markersize=3 ) 
        axes[i].fill_between(all_data[dp].q05['cases'][ti:tf].index, all_data[dp].q05['cases'][ti:tf], all_data[dp].q95['cases'][ti:tf], alpha = .5, color = 'darkblue', linewidth = 0)
        axes[i].plot(all_data[dp].q50['cases'][ti:tf], alpha = 1,linestyle='-', linewidth = 2, color = 'darkblue')
        axt.bar(pd.date_range(ti,tf, freq='W-SAT').date, rain[dp].resample('W-SAT').sum()[ti:tf], label = r'Rainfall', color = 'darkblue', width=7, alpha = 1)

        axes[i].set_title(dp)
        axes[i].set_ylim(0, 500)
        axt.set_ylim(1000, 0)
        #axt.set_ylim(4*max(rain[dp].resample('W-SAT').sum()[t_start:t_for]),0) # check if only reverse y
        axes[i].set_xlim(ti, tf)
        if i%5 == 4:
            axt.set_ylabel('Rainfall [mm/week]')
            axes[i].get_yaxis().set_visible(False)
        
        elif i%5 == 0:
            axes[i].set_ylabel('Reported cholera cases')
        if i%5 != 4:
            axt.get_yaxis().set_visible(False)

        if (dp not in scenario.not_dep) and (scenario_str != 'S0'):
            # convert to matplotlib date representation
            start = mdates.date2num(scenario.t_vacc_start[dp])
            end = mdates.date2num(scenario.t_vacc_end[dp])
            width = end - start
            rect = Rectangle((start, 0), width, 1000+max(all_data[dp].q95['cases']), color='orange', alpha= 0.1)
            axes[i].add_patch(rect) 
            axes[i].add_artist(rect)
            rx, ry = rect.get_xy()
            cx = rx + rect.get_width()/2.0
            cy = ry + rect.get_height()/1.5
    
    for ax in axes:
        ax.label_outer()

    fig.autofmt_xdate()
    fig.tight_layout()

    plt.savefig('output/Results/' + scenario_str + '.png', bbox_inches='tight')
    

def save_result(to_save, scenario_str, folder_name = ''):
    csv_all_q50 =  pd.DataFrame(0, index = pd.DatetimeIndex(start =  t_start, 
                                      end = t_for, freq = 'W-SAT'), columns = dept_name)

    csv_all_q05 =  pd.DataFrame(0, index = pd.DatetimeIndex(start =  t_start, 
                                      end = t_for, freq = 'W-SAT'), columns = dept_name)

    csv_all_q95 =  pd.DataFrame(0, index = pd.DatetimeIndex(start =  t_start, 
                                      end = t_for, freq = 'W-SAT'), columns = dept_name)


    for dp in departements:
        csv_all_q50[dp] = to_save[dp].q50['cases']
        csv_all_q95[dp] = to_save[dp].q95['cases']
        csv_all_q05[dp] = to_save[dp].q05['cases']
        
    csv_all_q05.to_csv(folder_name + scenario_str + '_q05.csv', index_label='date')
    csv_all_q95.to_csv(folder_name + scenario_str + '_q95.csv', index_label='date')
    csv_all_q50.to_csv(folder_name + scenario_str + '_q50.csv', index_label='date')
    pd.concat([csv_all_q05.sum(axis=1), csv_all_q50.sum(axis=1), csv_all_q95.sum(axis=1)], 
              axis = 1, keys = ['q05', 'q50', 'q95']).to_csv(folder_name + scenario_str + '_national.csv', index_label='date')

    for dp in departements:
        csv_all_q50[dp] = to_save[dp].q50['C']
        csv_all_q95[dp] = to_save[dp].q95['C']
        csv_all_q05[dp] = to_save[dp].q05['C']
        
    csv_all_q05.to_csv(folder_name + scenario_str + '_q05_s.csv', index_label='date')
    csv_all_q95.to_csv(folder_name + scenario_str + '_q95_s.csv', index_label='date')
    csv_all_q50.to_csv(folder_name + scenario_str + '_q50_s.csv', index_label='date')
    pd.concat([csv_all_q05.sum(axis=1), csv_all_q50.sum(axis=1), csv_all_q95.sum(axis=1)], 
              axis = 1, keys = ['q05', 'q50', 'q95']).to_csv(folder_name + scenario_str + '_national_s.csv', index_label='date')

scenarios_df = pd.read_csv('haiti-data/scenarios.csv', index_col = 0)

scenarios = {}
for sid, row in scenarios_df.iterrows():
    not_dep = []
    course_year = 2
    if (row['Roll-out'] == 2):
        not_dep = ['Ouest','Nord-Ouest','Sud', 'Nippes','Nord-Est', 'Sud-Est','Grande_Anse']
    elif (row['Roll-out'] == 3):
        course_year = 5
    elif (row['Roll-out'] == 4):
        not_dep = ['Nord-Ouest','Sud', 'Nippes','Nord-Est', 'Sud-Est','Grande_Anse']
        
    percent_completely_unvaccinated = 0
    percent_onedose = 0
    percent_twodoses = 0
    
    if (row['Coverage'] == 1):
        percent_completely_unvaccinated = 20
        percent_onedose = 10
        percent_twodoses = 70
    elif (row['Coverage'] == 2):
        percent_completely_unvaccinated = 40
        percent_onedose = 20
        percent_twodoses = 40
    elif (row['Coverage'] == 3):
        percent_completely_unvaccinated = 3.33
        percent_onedose = 1.67
        percent_twodoses = 95
    ve = row['VE']
            
    #if (row['Priority'] == 1):
    if True:
        scenarios['S' + str(sid)] = VaccinationScenario(course_year, 
                                                      percent_completely_unvaccinated, 
                                                      percent_onedose, 
                                                      percent_twodoses,
                                                      ve,
                                                      not_dep)


S0 = VaccinationScenario(50, 
                         99.9999999, 
                         0.00000001, 
                         0.00,
                         ve = 1)
scenarios['S0'] = S0


with mp.Pool(processes=n_proc) as pool:
    for _ in pool.imap_unordered(run_sim, scenarios):
        pass