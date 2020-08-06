import numpy as np
import json
import pandas as pd
import random
import matplotlib.pyplot as plt
import seaborn as sns

# Overal constant data
TS_list = [1,2,3,4]
RS_list = [1,2,3] # OPT / DEST / HDG
GF_list = [False, True, False, True]
wind_list = [False, False, True, True]

def overal_scenario_statistics_TS_RS(TS, RS, gf_defined = False):
    # Create dictionary to put statistics in
    statistics_dict = {}

    # Load scenario data
    scenario_report_json = open("thesis_tools/results/reports/scenario_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
    scenario_data = json.load(scenario_report_json)
    scenario_report_json.close()

    # Load conflict data
    conflict_report_json = open("thesis_tools/results/reports/conflict_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
    conflict_data = json.load(conflict_report_json)
    conflict_report_json.close()

    # initialize counters
    n_scenarios = len(scenario_data)
    n_conflicts = 0
    n_PZ_violated = 0
    n_PZ_violated_scenario = 0
    n_GF_violated = 0
    n_GF_violated_scenario = 0

    for i in range(n_scenarios):
        # filter conflict data for scenario number
        conflict_data_scenario = conflict_data[n_conflicts : n_conflicts + scenario_data[i]['n_conflicts']]

        n_conflicts += scenario_data[i]['n_conflicts'] 
        n_PZ_violated += scenario_data[i]['n_PZ_violated']
        if (scenario_data[i]['n_PZ_violated'] > 0):
            n_PZ_violated_scenario += 1
        if (gf_defined):
            n_GF_violated += scenario_data[i]['gf_violated0'] + scenario_data[i]['gf_violated1']
            if (scenario_data[i]['gf_violated0'] > 0):
                n_GF_violated_scenario += 1
            if (scenario_data[i]['gf_violated1'] > 0):
                n_GF_violated_scenario += 1

    t_conflicts = []
    for i in range(len(conflict_data)):
        t_conflicts.append(conflict_data[i]['t_conflict'])
    
    t_conflicts_avg = np.average(t_conflicts)

    statistics_dict['n_conflicts_avg'] =  n_conflicts / len(scenario_data)
    statistics_dict['t_conflicts_avg'] = t_conflicts_avg
    statistics_dict['C_IPR_avg'] = (n_conflicts - n_PZ_violated) / n_conflicts
    statistics_dict['S_IPR_avg'] = (n_scenarios - n_PZ_violated_scenario) / n_scenarios
    if (gf_defined):
        statistics_dict['C_VPRG_avg'] = (2. * n_conflicts - n_GF_violated) / (2. * n_conflicts)
        statistics_dict['S_VPRG_avg'] = (2. * n_scenarios - n_GF_violated_scenario) / (2. * n_scenarios)
    return statistics_dict

def create_scenario_statistics_df():
    # Create empty data dictionary
    data_dict = {   'TS':[], \
                    'RS':[], \
                    'n_conflicts_avg':[], \
                    't_conflicts_avg':[], \
                    'C_IPR':[], \
                    'S_IPR':[], \
                    'C_VPRG':[], \
                    'S_VPRG':[] }

    # Loop through Test Series and Rulesets
    idx = 0
    for TS in TS_list:
        for RS in RS_list:
            gf_defined = GF_list[idx]
            statistics = overal_scenario_statistics_TS_RS(TS, RS, gf_defined)

            # Add data to data dictionary
            data_dict['TS'].append(TS)
            data_dict['RS'].append(RS)
            data_dict['n_conflicts_avg'].append(statistics['n_conflicts_avg'])
            data_dict['t_conflicts_avg'].append(statistics['t_conflicts_avg'])            
            data_dict['C_IPR'].append(statistics['C_IPR_avg'])
            data_dict['S_IPR'].append(statistics['S_IPR_avg'])
            if gf_defined:
                data_dict['C_VPRG'].append(statistics['C_VPRG_avg'])
                data_dict['S_VPRG'].append(statistics['S_VPRG_avg'])
            else:
                data_dict['C_VPRG'].append(np.nan)
                data_dict['S_VPRG'].append(np.nan)
        idx += 1

    df = pd.DataFrame(data_dict)
    return df

#######################################################
# Hypothesis 1: Effect of implementation of geofences #
#######################################################

def hypothesis1_data_generator(sample_size_conflicts, sample_size_scenarios):
    data_IPR_C = []
    data_IPR_S = []
    data_VPRG_C = []
    data_VPRG_S = []

    idx = 0
    for TS in TS_list:
        GF_defined = GF_list[idx]
        for RS in RS_list:
            # Load scenario data
            scenario_report_json = open("thesis_tools/results/reports/scenario_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            scenario_data = json.load(scenario_report_json)
            scenario_report_json.close()

            # Load conflict data
            conflict_report_json = open("thesis_tools/results/reports/conflict_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            conflict_data = json.load(conflict_report_json)
            conflict_report_json.close()

            # Create randomized scenario data array
            scenario_data_randomised = scenario_data
            random.shuffle(scenario_data_randomised)

            # Create randomized conflict data array
            conflict_data_randomised = conflict_data
            random.shuffle(conflict_data_randomised)

            # Number of scenarios and conflicts
            n_scenarios = len(scenario_data)
            n_conflicts = len(conflict_data)

            # Analyse data on conflict base
            for i in range (int(n_conflicts / sample_size_conflicts)):
                n_intrusions = 0
                n_violations = 0
                for j in range(sample_size_conflicts):
                    idx_conflict = i * sample_size_conflicts + j
                    n_intrusions += conflict_data_randomised[idx_conflict]['intrusion']
                    if (GF_defined):
                        n_violations += conflict_data_randomised[idx_conflict]['violated_gf0']
                        n_violations += conflict_data_randomised[idx_conflict]['violated_gf1']
                IPR_C = (sample_size_conflicts - n_intrusions) / sample_size_conflicts
                data_IPR_C.append(['TS' + str(TS), 'RS' + str(RS), IPR_C])
                if (GF_defined):
                    VPRG_C = (2. * sample_size_conflicts - n_violations) / (2. * sample_size_conflicts)
                    data_VPRG_C.append(['TS' + str(TS), 'RS' + str(RS), VPRG_C])

            # Analyse data on scenario base
            for i in range(int(n_scenarios / sample_size_scenarios)):
                n_intrusions = 0 # Number of scenarios have intrusions
                n_violations = 0 # Number of UAVs involved in scenario have violated the GF
                for j in range(sample_size_scenarios):
                    idx_scenario = i * sample_size_scenarios + j
                    n_intrusions += scenario_data_randomised[idx_scenario]['n_PZ_violated'] > 0
                    if (GF_defined):
                        n_violations += scenario_data_randomised[idx_scenario]['gf_violated0'] > 0
                        n_violations += scenario_data_randomised[idx_scenario]['gf_violated1'] > 0
                IPR_S = (sample_size_scenarios - n_intrusions) / sample_size_scenarios
                data_IPR_S.append(['TS' + str(TS), 'RS' + str(RS), IPR_S])
                if (GF_defined):
                    VPRG_S = (2. * sample_size_scenarios - n_violations) / (2. * sample_size_scenarios)
                    data_VPRG_S.append(['TS' + str(TS), 'RS' + str(RS), VPRG_S])
        idx += 1
    
    # write dataframes to files
    df_IPR_C = pd.DataFrame(data_IPR_C, columns = ['TS', 'RS', 'IPR_C'])
    df_IPR_S = pd.DataFrame(data_IPR_S, columns = ['TS', 'RS', 'IPR_S'])
    df_VPRG_C = pd.DataFrame(data_VPRG_C, columns = ['TS', 'RS', 'VPRG_C'])
    df_VPRG_S = pd.DataFrame(data_VPRG_S, columns = ['TS', 'RS', 'VPRG_S'])

    df_IPR_C.to_csv("thesis_tools/results/performance/hypothesis1_ipr_c.csv")
    df_IPR_S.to_csv("thesis_tools/results/performance/hypothesis1_ipr_s.csv")
    df_VPRG_C.to_csv("thesis_tools/results/performance/hypothesis1_vprg_c.csv")
    df_VPRG_S.to_csv("thesis_tools/results/performance/hypothesis1_vprg_s.csv")

def hypothesis1_IPR_box_whiskerplot_creator():
    df_IPR_C = pd.read_csv("thesis_tools/results/performance/hypothesis1_ipr_c.csv")
    df_IPR_S = pd.read_csv("thesis_tools/results/performance/hypothesis1_ipr_s.csv")

    df_IPR_C['TS'] = df_IPR_C['TS'].map({'TS1': 'no geofence\nno wind', 'TS2': 'geofence\nno wind', 'TS3': 'no geofence\nwind', 'TS4': 'geofence\nwind'})
    df_IPR_C['RS'] = df_IPR_C['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})

    df_IPR_S['TS'] = df_IPR_S['TS'].map({'TS1': 'no geofence\nno wind', 'TS2': 'geofence\nno wind', 'TS3': 'no geofence\nwind', 'TS4': 'geofence\nwind'})
    df_IPR_S['RS'] = df_IPR_S['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})

    plt.figure()
    sns.boxplot(x='TS',y='IPR_C', hue="RS", data=df_IPR_C, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.figure()
    sns.boxplot(x='TS',y='IPR_S', hue="RS", data=df_IPR_S, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.show()


def hypothesis1_VPRG_box_whiskerplot_creator():
    df_VPRG_C = pd.read_csv("thesis_tools/results/performance/hypothesis1_vprg_c.csv")
    df_VPRG_S = pd.read_csv("thesis_tools/results/performance/hypothesis1_vprg_s.csv")
    
    df_VPRG_C['TS'] = df_VPRG_C['TS'].map({'TS2': 'no wind', 'TS4': 'wind'})
    df_VPRG_C['RS'] = df_VPRG_C['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})
    
    df_VPRG_S['TS'] = df_VPRG_S['TS'].map({'TS2': 'no wind', 'TS4': 'wind'})
    df_VPRG_S['RS'] = df_VPRG_S['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})
    
    plt.figure()
    sns.boxplot(x='TS',y='VPRG_C', hue="RS", data=df_VPRG_C, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.figure()
    sns.boxplot(x='TS',y='VPRG_S', hue="RS", data=df_VPRG_S, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.show()

################################
# Hypothesis 2: Effect of wind #
################################

""" Only test series (3 and 4) with wind will be plotted for this hypothesis. 
    Wind will be subdivided in 3 categories (weak, medium, strong)
    So 4 plots (with, without geofence) for VPRG_C, VPRG_S, IPR_C, IPR_S"""

def hypothesis2_data_generator(sample_size_conflicts, sample_size_scenarios, first_limit, second_limit):
    # Create empty data arrays
    data_IPR_C = []
    data_IPR_S = []
    data_VPRG_C = []
    data_VPRG_S = []

    TS_indices = [2,3]
    idx = 0
    for TS in np.array(TS_list)[TS_indices]:
        GF_defined = np.array(GF_list)[TS_indices[idx]]
        for RS in RS_list:
            # Load scenario data
            scenario_report_json = open("thesis_tools/results/reports/scenario_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            scenario_data = json.load(scenario_report_json)
            scenario_report_json.close()

            # Load conflict data
            conflict_report_json = open("thesis_tools/results/reports/conflict_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            conflict_data = json.load(conflict_report_json)
            conflict_report_json.close()

            # Create randomized scenario data array
            scenario_data_randomised = scenario_data
            random.shuffle(scenario_data_randomised)

            # Create randomized conflict data array
            conflict_data_randomised = conflict_data
            random.shuffle(conflict_data_randomised)

            # Number of scenarios and conflicts
            n_scenarios = len(scenario_data)
            n_conflicts = len(conflict_data)

            # categorise scenario data by windspeed
            scenario_data_categorised = [[], [], []]
            
            for i in range(n_scenarios):
                if scenario_data_randomised[i]['windspeed'] < first_limit:
                    scenario_data_categorised[0].append(scenario_data_randomised[i])
                elif scenario_data_randomised[i]['windspeed'] < second_limit:
                    scenario_data_categorised[1].append(scenario_data_randomised[i])
                else:
                    scenario_data_categorised[2].append(scenario_data_randomised[i])

            # categorise conflict data by windspeed
            conflict_data_categorised = [[], [], []]

            for i in range(n_conflicts):
                if conflict_data_randomised[i]['windspeed'] < first_limit:
                    conflict_data_categorised[0].append(conflict_data_randomised[i])
                elif conflict_data_randomised[i]['windspeed'] < second_limit:
                    conflict_data_categorised[1].append(conflict_data_randomised[i])
                else:
                    conflict_data_categorised[2].append(conflict_data_randomised[i])

            # Analyse data on conflict base:
            for i in range(3):
                n_conflicts_category = len(conflict_data_categorised[i])
                wind_category = ''
                if i == 0: 
                    wind_category = 'low'
                elif i == 1:
                    wind_category = 'medium'
                else:
                    wind_category = 'strong'

                for j in range(int(n_conflicts_category / sample_size_conflicts)):
                    n_intrusions = 0
                    n_violations = 0
                    for k in range(sample_size_conflicts):
                        idx_conflict = j * sample_size_conflicts + k
                        n_intrusions += conflict_data_categorised[i][idx_conflict]['intrusion']
                        if GF_defined:
                            n_violations += conflict_data_categorised[i][idx_conflict]['violated_gf0']
                            n_violations += conflict_data_categorised[i][idx_conflict]['violated_gf1']
                    IPR_C = (sample_size_conflicts - n_intrusions) / sample_size_conflicts
                    data_IPR_C.append(['TS' + str(TS), 'RS' + str(RS), IPR_C, wind_category])
                    if (GF_defined):
                        VPRG_C = (2. * sample_size_conflicts - n_violations) / (2. * sample_size_conflicts)
                        data_VPRG_C.append(['TS' + str(TS), 'RS' + str(RS), VPRG_C, wind_category])

            # Analyse data on scenario base:
            for i in range(3):
                n_scenarios_category = len(scenario_data_categorised[i])
                wind_category = ''
                if i == 0: 
                    wind_category = 'low'
                elif i == 1:
                    wind_category = 'medium'
                else:
                    wind_category = 'strong'

                for j in range(int(n_scenarios_category / sample_size_scenarios)):
                    n_intrusions = 0
                    n_violations = 0
                    for k in range(sample_size_scenarios):
                        idx_scenario = j * sample_size_scenarios + k
                        n_intrusions += scenario_data_categorised[i][idx_scenario]['n_PZ_violated'] > 0
                        if GF_defined:
                            n_violations += scenario_data_categorised[i][idx_scenario]['gf_violated0'] > 0
                            n_violations += scenario_data_categorised[i][idx_scenario]['gf_violated1'] > 0
                    IPR_S = (sample_size_scenarios - n_intrusions) / sample_size_scenarios
                    data_IPR_S.append(['TS' + str(TS), 'RS' + str(RS), IPR_S, wind_category])
                    if (GF_defined):
                        VPRG_S = (2. * sample_size_scenarios - n_violations) / (2. * sample_size_scenarios)
                        data_VPRG_S.append(['TS' + str(TS), 'RS' + str(RS), VPRG_S, wind_category])
        idx += 1

    # write dataframes to files
    df_IPR_C = pd.DataFrame(data_IPR_C, columns = ['TS', 'RS', 'IPR_C', 'wind'])
    df_IPR_S = pd.DataFrame(data_IPR_S, columns = ['TS', 'RS', 'IPR_S', 'wind'])
    df_VPRG_C = pd.DataFrame(data_VPRG_C, columns = ['TS', 'RS', 'VPRG_C', 'wind'])
    df_VPRG_S = pd.DataFrame(data_VPRG_S, columns = ['TS', 'RS', 'VPRG_S', 'wind'])

    df_IPR_C.to_csv("thesis_tools/results/performance/hypothesis2_ipr_c.csv")
    df_IPR_S.to_csv("thesis_tools/results/performance/hypothesis2_ipr_s.csv")
    df_VPRG_C.to_csv("thesis_tools/results/performance/hypothesis2_vprg_c.csv")
    df_VPRG_S.to_csv("thesis_tools/results/performance/hypothesis2_vprg_s.csv")

def hypothesis2_IPR_box_whiskerpplot_creator():
    df_IPR_C = pd.read_csv("thesis_tools/results/performance/hypothesis2_ipr_c.csv")
    df_IPR_S = pd.read_csv("thesis_tools/results/performance/hypothesis2_ipr_s.csv")

    df_IPR_C['RS'] = df_IPR_C['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})
    df_IPR_S['RS'] = df_IPR_S['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})

    df_IPR_C_nogf = df_IPR_C[df_IPR_C['TS'] == 'TS3']
    df_IPR_S_nogf = df_IPR_S[df_IPR_S['TS'] == 'TS3']
    df_IPR_C_gf = df_IPR_C[df_IPR_C['TS'] == 'TS4']
    df_IPR_S_gf = df_IPR_S[df_IPR_S['TS'] == 'TS4']

    plt.figure()
    sns.boxplot(x='RS',y='IPR_C', hue='wind', data=df_IPR_C_nogf, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.figure()
    sns.boxplot(x='RS',y='IPR_C', hue='wind', data=df_IPR_C_gf, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.figure()
    sns.boxplot(x='RS',y='IPR_S', hue='wind', data=df_IPR_S_nogf, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.figure()
    sns.boxplot(x='RS',y='IPR_S', hue='wind', data=df_IPR_S_gf, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.show()

def hypothesis2_VPRG_box_whiskerpplot_creator():
    df_VPRG_C = pd.read_csv("thesis_tools/results/performance/hypothesis2_vprg_c.csv")
    df_VPRG_S = pd.read_csv("thesis_tools/results/performance/hypothesis2_vprg_s.csv")

    df_VPRG_C['RS'] = df_VPRG_C['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})
    df_VPRG_S['RS'] = df_VPRG_S['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})

    plt.figure()
    sns.boxplot(x='RS',y='VPRG_C', hue='wind', data=df_VPRG_C, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.figure()
    sns.boxplot(x='RS',y='VPRG_S', hue='wind', data=df_VPRG_S, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.show()


#################################################
# Hypothesis 3: Effect of distance wrt geofence #
#################################################

""" Only test series (2 and 4) with geofences will be plotted for this hypothesis. 
    Distance to geofence will be subdivided in 3 categories (weak, medium, strong)
    So 4 plots (without geofence) for VPRG_C, VPRG_S, IPR_C, IPR_S"""

def hypothesis3_data_generator(sample_size_conflicts, sample_size_scenarios, first_limit, second_limit):
    # Create empty data arrays
    data_IPR_C = []
    data_IPR_S = []
    data_VPRG_C = []
    data_VPRG_S = []

    TS_indices = [1,3]
    for TS in np.array(TS_list)[TS_indices]:
        for RS in RS_list:
            # Load scenario data
            scenario_report_json = open("thesis_tools/results/reports/scenario_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            scenario_data = json.load(scenario_report_json)
            scenario_report_json.close()

            # Load conflict data
            conflict_report_json = open("thesis_tools/results/reports/conflict_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            conflict_data = json.load(conflict_report_json)
            conflict_report_json.close()

            # Number of scenarios and conflicts
            n_scenarios = len(scenario_data)
            n_conflicts = len(conflict_data)

            # split conflict data dependent on dist to gf per UAV involved in each conflict
            conflict_data_splitted = []

            for i in range(n_conflicts):
                conflict_gf0 = {'intrusion': conflict_data[i]['intrusion'],
                                'violation': conflict_data[i]['violated_gf0'],
                                'dist_gf': conflict_data[i]['dist_gf0']}

                conflict_gf1 = {'intrusion': conflict_data[i]['intrusion'],
                                'violation': conflict_data[i]['violated_gf1'],
                                'dist_gf': conflict_data[i]['dist_gf1']}

                conflict_data_splitted.append(conflict_gf0)
                conflict_data_splitted.append(conflict_gf1)

            # randomize order of splitted conflict data
            conflict_data_randomised = conflict_data_splitted
            random.shuffle(conflict_data_randomised)

            # split scenario data dependen on dist to gf per UAV involved in each conflict
            scenario_data_splitted = []

            for i in range(n_scenarios):
                scenario_gf0 = {'intrusion': scenario_data[i]['n_PZ_violated'] > 0,
                                'violation': scenario_data[i]['gf_violated0'],
                                'dist_gf': scenario_data[i]['min_pos_dist_gf0']}

                scenario_gf1 = {'intrusion': scenario_data[i]['n_PZ_violated'] > 0,
                                'violation': scenario_data[i]['gf_violated1'],
                                'dist_gf': scenario_data[i]['min_pos_dist_gf1']}

                scenario_data_splitted.append(scenario_gf0)
                scenario_data_splitted.append(scenario_gf1)

            # randomize order of splitted scenario data
            scenario_data_randomised = scenario_data_splitted
            random.shuffle(scenario_data_randomised)

            # categorise conflict data by dist to geofence
            conflict_data_categorised = [[], [], []]

            for i in range(2 * n_conflicts):
                # first filer out conflicts that are outside geofence:
                if conflict_data_randomised[i]['dist_gf'] < 0:
                    pass

                if conflict_data_randomised[i]['dist_gf'] < first_limit:
                    conflict_data_categorised[0].append(conflict_data_randomised[i])
                elif conflict_data_randomised[i]['dist_gf'] < second_limit:
                    conflict_data_categorised[1].append(conflict_data_randomised[i])
                else:
                    conflict_data_categorised[2].append(conflict_data_randomised[i])

            # categorise scenario data by dist to geofence
            scenario_data_categorised = [[], [], []]

            for i in range(2 * n_scenarios):
                if scenario_data_randomised[i]['dist_gf'] < first_limit:
                    scenario_data_categorised[0].append(scenario_data_randomised[i])
                elif scenario_data_randomised[i]['dist_gf'] < second_limit:
                    scenario_data_categorised[1].append(scenario_data_randomised[i])
                else:
                    scenario_data_categorised[2].append(scenario_data_randomised[i])

            # Analyse data on conflict base:
            for i in range(3):
                n_conflicts_category = len(conflict_data_categorised[i])
                dist_category = ''
                if i == 0:
                    dist_category = 'small'
                elif i == 1:
                    dist_category = 'medium'
                else:
                    dist_category = 'large'
                
                for j in range(int(n_conflicts_category / sample_size_conflicts)):
                    n_intrusions = 0
                    n_violations = 0
                    for k in range(sample_size_conflicts):
                        idx_conflict = j * sample_size_conflicts + k
                        n_intrusions += conflict_data_categorised[i][idx_conflict]['intrusion']
                        n_violations += conflict_data_categorised[i][idx_conflict]['violation']
                    IPR_C = (sample_size_conflicts - n_intrusions) / sample_size_conflicts
                    data_IPR_C.append(['TS' + str(TS), 'RS' + str(RS), IPR_C, dist_category])
                    VPRG_C = (sample_size_conflicts - n_violations) / (sample_size_conflicts)
                    data_VPRG_C.append(['TS' + str(TS), 'RS' + str(RS), VPRG_C, dist_category])
            
            # Analyse data on scenario base:
            for i in range(3):
                n_scenarios_category = len(scenario_data_categorised[i])
                dist_category = ''
                if i == 0:
                    dist_category = 'small'
                elif i == 1:
                    dist_category = 'medium'
                else:
                    dist_category = 'large'
                
                for j in range(int(n_scenarios_category / sample_size_scenarios)):
                    n_intrusions = 0
                    n_violations = 0
                    for k in range(sample_size_scenarios):
                        idx_scenario = j * sample_size_scenarios + k
                        n_intrusions += scenario_data_categorised[i][idx_scenario]['intrusion']
                        n_violations += scenario_data_categorised[i][idx_scenario]['violation']
                    IPR_S = (sample_size_scenarios - n_intrusions) / sample_size_scenarios
                    data_IPR_S.append(['TS' + str(TS), 'RS' + str(RS), IPR_S, dist_category])
                    VPRG_S = (sample_size_scenarios - n_violations) / (sample_size_scenarios)
                    data_VPRG_S.append(['TS' + str(TS), 'RS' + str(RS), VPRG_S, dist_category])
    
    # write dataframes to files
    df_IPR_C = pd.DataFrame(data_IPR_C, columns = ['TS', 'RS', 'IPR_C', 'dist'])
    df_IPR_S = pd.DataFrame(data_IPR_S, columns = ['TS', 'RS', 'IPR_S', 'dist'])
    df_VPRG_C = pd.DataFrame(data_VPRG_C, columns = ['TS', 'RS', 'VPRG_C', 'dist'])
    df_VPRG_S = pd.DataFrame(data_VPRG_S, columns = ['TS', 'RS', 'VPRG_S', 'dist'])

    df_IPR_C.to_csv("thesis_tools/results/performance/hypothesis3_ipr_c.csv")
    df_IPR_S.to_csv("thesis_tools/results/performance/hypothesis3_ipr_s.csv")
    df_VPRG_C.to_csv("thesis_tools/results/performance/hypothesis3_vprg_c.csv")
    df_VPRG_S.to_csv("thesis_tools/results/performance/hypothesis3_vprg_s.csv")

def hypothesis3_IPR_box_whiskerplot_creator():
    df_IPR_C = pd.read_csv("thesis_tools/results/performance/hypothesis3_ipr_c.csv")
    df_IPR_S = pd.read_csv("thesis_tools/results/performance/hypothesis3_ipr_s.csv")

    df_IPR_C['RS'] = df_IPR_C['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})
    df_IPR_S['RS'] = df_IPR_S['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})

    df_IPR_C_nowind = df_IPR_C[df_IPR_C['TS'] == 'TS2']
    df_IPR_S_nowind = df_IPR_S[df_IPR_S['TS'] == 'TS2']
    df_IPR_C_wind = df_IPR_C[df_IPR_C['TS'] == 'TS4']
    df_IPR_S_wind = df_IPR_S[df_IPR_S['TS'] == 'TS4']

    plt.figure()
    sns.boxplot(x='RS',y='IPR_C', hue='dist', data=df_IPR_C_nowind, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.figure()
    sns.boxplot(x='RS',y='IPR_C', hue='dist', data=df_IPR_C_wind, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.figure()
    sns.boxplot(x='RS',y='IPR_S', hue='dist', data=df_IPR_S_nowind, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.figure()
    sns.boxplot(x='RS',y='IPR_S', hue='dist', data=df_IPR_S_wind, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.show()

def hypothesis3_VPRG_box_whiskerplot_creator():
    df_VPRG_C = pd.read_csv("thesis_tools/results/performance/hypothesis3_vprg_c.csv")
    df_VPRG_S = pd.read_csv("thesis_tools/results/performance/hypothesis3_vprg_S.csv")

    df_VPRG_C['RS'] = df_VPRG_C['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})
    df_VPRG_S['RS'] = df_VPRG_S['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})

    df_VPRG_C_nowind = df_VPRG_C[df_VPRG_C['TS'] == 'TS2']
    df_VPRG_S_nowind = df_VPRG_S[df_VPRG_S['TS'] == 'TS2']
    df_VPRG_C_wind = df_VPRG_C[df_VPRG_C['TS'] == 'TS4']
    df_VPRG_S_wind = df_VPRG_S[df_VPRG_S['TS'] == 'TS4']

    plt.figure()
    sns.boxplot(x='RS',y='VPRG_C', hue='dist', data=df_VPRG_C_nowind, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.figure()
    sns.boxplot(x='RS',y='VPRG_C', hue='dist', data=df_VPRG_C_wind, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.figure()
    sns.boxplot(x='RS',y='VPRG_S', hue='dist', data=df_VPRG_S_nowind, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.figure()
    sns.boxplot(x='RS',y='VPRG_S', hue='dist', data=df_VPRG_S_wind, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.6)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.show()

#print(create_scenario_statistics_df())
#hypothesis1_data_generator(200, 200)
#hypothesis1_IPR_box_whiskerplot_creator()
#hypothesis1_VPRG_box_whiskerplot_creator()
#hypothesis2_data_generator(200,200,5,15)
#hypothesis2_IPR_box_whiskerpplot_creator()
#hypothesis2_VPRG_box_whiskerpplot_creator()
#hypothesis3_data_generator(200,200,100,150)
#hypothesis3_IPR_box_whiskerplot_creator()
#hypothesis3_VPRG_box_whiskerplot_creator()