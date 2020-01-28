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

output_dir = 'output_16-04-init/'
nsim = 500
run_lvl = 3
n_proc = 10 
        
dept_avail = os.listdir(output_dir)

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
rainfall = pd.read_csv('haiti-data/fromAzman/rainfall.csv', index_col = 0, parse_dates = True)
cases    = pd.read_csv('haiti-data/fromAzman/cases_corrected.csv', index_col=0, parse_dates =True)
rain     = pd.read_csv('haiti-data/proj/rainfall.csv', index_col = 0, parse_dates = True)

compartments = ['S', 'I', 'A', 'RA1', 'RA2', 'RA3', 'RI1', 'RI2', 'RI3', 'W', 'B', 'cases', 'C',
                 "VSd", "VRI1d", "VRI2d", "VRI3d", "VRA1d", "VRA2d", "VRA3d",
                 "VSdd", "VRI1dd", "VRI2dd", "VRI3dd", "VRA1dd", "VRA2dd", "VRA3dd",
                 "VSd_alt", "VRI1d_alt", "VRI2d_alt", "VRI3d_alt", "VRA1d_alt", "VRA2d_alt", "VRA3d_alt",
                 "VSdd_alt", "VRI1dd_alt", "VRI2dd_alt", "VRI3dd_alt", "VRA1dd_alt", "VRA2dd_alt", "VRA3dd_alt"]

stream = open('haiti-data/input_parameters.yaml', 'r')
input_parameters = yaml.load(stream)

dept_name = [list(pop.keys())[0] for pop in input_parameters['population']]

t_start = input_parameters['t_start']
t_for = datetime.date(2029,12,20)

index = pd.DatetimeIndex(start =  t_start,  end = t_for, freq = 'W-SAT')

def make_genuine_mobility_file():
    covar_init = pd.concat([cases[t_start:]]*6, ignore_index=True)[0:-506]
    covar_init.index = pd.DatetimeIndex(start =  datetime.date(2010,10,23), 
                                      end = t_for, freq = 'W-SAT')
    covar_init.to_csv('covar_mob.csv', index_label='date')
    return covar_init

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
def run_sim(dp):
    # see second answer of https://stackoverflow.com/questions/25175530/can-rpy2-code-be-run-in-parallel on why it starts different 
    # R instances
    r_source = robjects.r['source'];
    dept_data = {}
    r_options = robjects.r['options']
    r_options(warn=-1)
    robjects.r('departement  <- "' + dp + '"')
    robjects.r('output_dir   <- "' + output_dir + '"')
    robjects.r('scenario     <- "' + scenario_str + '"')
    robjects.r('run_level    <- '  + str(run_lvl))
    robjects.r('nsim         <- '  + str(nsim))
    robjects.r('t_vacc_start <- "' + str(scenario.t_vacc_start[dp]) + '"')
    robjects.r('t_vacc_end   <- "' + str(scenario.t_vacc_end[dp]) + '"')
    robjects.r('p1d_reg      <- '  + str(scenario.p1d_reg[dp]))
    robjects.r('r_v_year     <- '  + str(scenario.r_v_year[dp]))
    robjects.r('cases_ext    <- '  + str(scenario.ve))
    r_source('./scripts/forecast_haitiOCV_mob.R')
    for comp in compartments:
        #temp = pandas2ri.ri2py(robjects.r[comp])
        temp = robjects.r[comp]
        temp.index = index
        temp.drop('date',axis=1, inplace = True)
        dept_data[comp] = temp
    return (dp, dept_data, robjects.r['params'], robjects.r('names(params)'))

def save_result(to_save, scenario_str, folder_name = ''):
    csv_all_q50 =  pd.DataFrame(0, index = pd.DatetimeIndex(start =  t_start, 
                                      end = t_for, freq = 'W-SAT'), columns = dept_name)

    csv_all_q05 =  pd.DataFrame(0, index = pd.DatetimeIndex(start =  t_start, 
                                      end = t_for, freq = 'W-SAT'), columns = dept_name)

    csv_all_q95 =  pd.DataFrame(0, index = pd.DatetimeIndex(start =  t_start, 
                                      end = t_for, freq = 'W-SAT'), columns = dept_name)


    for dp in dept_avail:
        csv_all_q50[dp] = to_save[dp]['cases']['q50']
        csv_all_q95[dp] = to_save[dp]['cases']['q95']
        csv_all_q05[dp] = to_save[dp]['cases']['q05']
        
    csv_all_q05.to_csv(folder_name + scenario_str + '_q05.csv', index_label='date')
    csv_all_q95.to_csv(folder_name + scenario_str + '_q95.csv', index_label='date')
    csv_all_q50.to_csv(folder_name + scenario_str + '_q50.csv', index_label='date')
    pd.concat([csv_all_q05.sum(axis=1), csv_all_q50.sum(axis=1), csv_all_q95.sum(axis=1)], 
              axis = 1, keys = ['q05', 'q50', 'q95']).to_csv(folder_name + scenario_str + '_national.csv', index_label='date')

    for dp in dept_avail:
        csv_all_q50[dp] = to_save[dp]['C']['q50']
        csv_all_q95[dp] = to_save[dp]['C']['q95']
        csv_all_q05[dp] = to_save[dp]['C']['q05']
        
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



for scenario_str, scenario  in scenarios.items():
    print ("Working on scenario " + scenario_str)
    
    # Projection without mobility
    all_data_vacc = {}
    make_genuine_mobility_file()
    with mp.Pool(processes=n_proc) as pool:
        for dp, dept_data, _, _ in pool.imap_unordered(run_sim, dept_avail):
            all_data_vacc[dp] = dept_data
            
    # Make mobility file
    vacc_init = pd.DataFrame(np.nan, pd.DatetimeIndex(start =  datetime.date(2010,10,23), 
                                      end = t_for, freq = 'W-SAT'), columns = dept_name)
    cases_df =  pd.DataFrame(np.nan, index = pd.DatetimeIndex(start =  t_start, 
                                      end = t_for, freq = 'W-SAT'), columns = dept_name)
    covar_mob_novacc = make_genuine_mobility_file()
    for dp in dept_avail:
        cases_df[dp] = all_data_vacc[dp]['cases']['q50'][scenario.t_vacc_start[dp]+datetime.timedelta(days = 90) : scenario.t_vacc_start[dp] + datetime.timedelta(days = 200)]
        if (sum(all_data_vacc[dp]['cases']['q50'][scenario.t_vacc_start[dp]: scenario.t_vacc_start[dp] + datetime.timedelta(days = 200)] == 0) >= 3):
        # No more cases at some point
            cases_df[dp] = 0
    for dp in dept_name:
        if dp in dept_avail:
            vacc_init[dp] = covar_mob_novacc[dp][datetime.date(2010,10,23):scenario.t_vacc_start[dp]]   # 
            vacc_init[dp].fillna(cases_df[dp].mean(), inplace = True)
        else:
            vacc_init[dp] = covar_mob_novacc[dp]  #0
    vacc_init.to_csv('vacc_init' + scenario_str +  '.csv', index_label='date')
    
    # Run projection with mobility
    vacc_init = pd.read_csv('vacc_init' + scenario_str +  '.csv', index_col = 0, parse_dates = True)
    vacc_init.to_csv('covar_mob.csv', index_label='date')
    
    all_data_vacc_mob = {}
    with mp.Pool(processes=n_proc) as pool:
        for dp, dept_data, _, _ in pool.imap_unordered(run_sim, dept_avail):
            all_data_vacc_mob[dp] = dept_data

    all_data_vacc_mob['Ouest']['cases'][:datetime.date(2015,7,1)] = np.nan
    
    # Save results
    dir_name = 'output/Results/' + scenario_str + '/'
    try:
        os.makedirs(dir_name)    
        print("Directory " , dir_name ,  " Created ")
    except FileExistsError:
        print("Directory " , dir_name ,  " already exists, writing into.")
    save_result(all_data_vacc_mob, scenario_str, folder_name = dir_name)
    
    # Do figure
    ti = input_parameters['t_start']
    tf = t_for

    fig, axes = plt.subplots((len(all_data_vacc_mob))//2, 2, figsize=(15,15), squeeze = True, dpi = 200)
    axes = axes.flatten()
    fig.patch.set_facecolor('white')

    for i, dp in enumerate(dept_avail):
        axt =  axes[i].twinx()
        axes[i].plot(cases[dp][ti:tf], marker='.', linestyle='-',color='k', linewidth=0, markersize=1 ) 
        axes[i].fill_between(all_data_vacc_mob[dp]['cases']['q05'][ti:tf].index, all_data_vacc_mob[dp]['cases']['q05'][ti:tf], all_data_vacc_mob[dp]['cases']['q95'][ti:tf], alpha = .5, color = 'red', linewidth = 0)
        axes[i].plot(all_data_vacc_mob[dp]['cases']['q50'][ti:tf], alpha = 1,linestyle='-', linewidth = 1, color = 'darkblue')
        axt.bar(pd.date_range(t_start,t_for, freq='W-SAT').date, rain[dp].resample('W-SAT').sum()[t_start:t_for], 
                label = r'Rainfall', color = 'darkblue', width=7, alpha = 1)
    
        axes[i].set_title(dp)
        axes[i].set_ylim(0)
        axt.set_ylim(2*max(rain[dp].resample('W-SAT').sum()[t_start:t_for]),0) # check if only reverse y
        axes[i].set_xlim(ti, tf)
        if dp not in scenario.not_dep:
            # convert to matplotlib date representation
            start = mdates.date2num(scenario.t_vacc_start[dp])
            end = mdates.date2num(scenario.t_vacc_end[dp])
            width = end - start
            rect = Rectangle((start, 0), width, 1000+max(all_data_vacc_mob[dp]['cases']['q95']), color='orange', alpha= 0.1)
            axes[i].add_patch(rect) 
            axes[i].add_artist(rect)
            rx, ry = rect.get_xy()
            cx = rx + rect.get_width()/2.0
            cy = ry + rect.get_height()/1.5

    fig.tight_layout()
    plt.savefig('output/Results/' + scenario_str + '.png', bbox_inches='tight')
