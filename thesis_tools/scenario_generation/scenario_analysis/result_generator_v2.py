import numpy as np
import json
import pandas as pd

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

print(create_scenario_statistics_df())