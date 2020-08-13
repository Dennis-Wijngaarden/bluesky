import numpy as np
import json
import random
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import scipy as sp

#Init randomizer
random.seed()

# Some parameters
sample_size = 200

# calculate IPR for each sample of a TS and RS
def IPR_sample_calculator(TS, RS, sample_size):
    IPRs = []

    conflict_report_json = open("thesis_tools/results/reports/conflict_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
    conflict_data = json.load(conflict_report_json)
    conflict_report_json.close()

    random.shuffle(conflict_data)

    # split in samples
    for i in range(int(len(conflict_data) / sample_size)):
        n_intrusion = 0 # Number of intrusions
        for j in range(sample_size):
            idx = i * sample_size + j
            n_intrusion += conflict_data[idx]['intrusion']
        IPR = (sample_size - n_intrusion) / sample_size
        IPRs.append(IPR)

    return IPRs

# Calculate SVPRG for each sample of a TS and RS (Scenario Violation Prevention Rate of the Geofence)
def SVPRG_sample_calculator(TS, RS, sample_size):
    SVPRGs = []

    scenario_report_json = open("thesis_tools/results/reports/scenario_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
    scenario_data = json.load(scenario_report_json)
    scenario_report_json.close()

    random.shuffle(scenario_data)

    # split om samples
    for i in range(int(len(scenario_data) / sample_size)):
        n_violations = 0
        for j in range(sample_size):
            idx = i * sample_size + j
            n_violations += scenario_data[idx]['gf_violated0']
            n_violations += scenario_data[idx]['gf_violated1']
        SVPRG = (2. * sample_size - n_violations) / (2. * sample_size)
        SVPRGs.append(SVPRG)
    
    return SVPRGs

# Calculate CVPRG for each sample of a TS and RS (Conflict Violation Prevention Rate of the Geofence)
def CVPRG_sample_calculator(TS, RS, sample_size):
    CVPRGs = []

    conflict_report_json = open("thesis_tools/results/reports/conflict_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
    conflict_data = json.load(conflict_report_json)
    conflict_report_json.close()

    random.shuffle(conflict_data)

    # split om samples
    for i in range(int(len(conflict_data) / sample_size)):
        n_violations = 0
        for j in range(sample_size):
            idx = i * sample_size + j
            n_violations += conflict_data[idx]['violated_gf0']
            n_violations += conflict_data[idx]['violated_gf1']
        CVPRG = (2. * sample_size - n_violations) / (2. * sample_size)
        CVPRGs.append(CVPRG)
    
    return CVPRGs

#############################
# PLOTTERS FOR HYPOTHESIS 1 #
#############################

def IPR_box_whiskerplot_creator():
    df = pd.read_csv("thesis_tools/results/performance/hypothesis1_ipr.csv")
    df['TS'] = df['TS'].map({'TS1': 'no geofence\nno wind', 'TS2': 'geofence\nno wind', 'TS3': 'no geofence\nwind', 'TS4': 'geofence\nwind'})
    df['RS'] = df['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})
    sns.boxplot(x='TS',y='IPR', hue="RS", data=df, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.85)
    plt.legend(ncol=3)
    plt.subplots_adjust(bottom=0.15)
    plt.show()
    return

def CVPRG_box_whiskerplot_creator():
    df = pd.read_csv("thesis_tools/results/performance/hypothesis1_vprg.csv")
    df['TS'] = df['TS'].map({'TS2': 'no wind', 'TS4': 'wind'})
    df['RS'] = df['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})
    sns.boxplot(x='TS',y='VPRG', hue="RS", data=df, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.85)
    plt.legend(ncol=3)
    plt.show()
    return

#############################
# PLOTTERS FOR HYPOTHESIS 2 #
#############################

def determine_wind_boundaries_for_plot():
    # check windspeeds in conflict files
    windspeeds = []
    for TS in [3, 4]:
        for RS in np.arange(1,4):
            conflict_report_json = open("thesis_tools/results/reports/conflict_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            conflict_data = json.load(conflict_report_json)
            conflict_report_json.close()
            for i in range(len(conflict_data)):
                windspeeds.append(conflict_data[i]['windspeed'])
    windspeeds.sort()
    first_boundary_idx = int(len(windspeeds)/3) 
    second_boundary_idx = int(len(windspeeds)//3 * 2)
    first_boundary = windspeeds[first_boundary_idx]
    second_boundary_idx = windspeeds[second_boundary_idx]
    max_windspeed = max(windspeeds)
    return first_boundary, second_boundary_idx, max_windspeed

def IPR_VPRG_wind_box_whiskerplot_creator():
    df_IPR = pd.read_csv("thesis_tools/results/performance/hypothesis2_ipr.csv")
    df_VPRG = pd.read_csv("thesis_tools/results/performance/hypothesis2_vprg.csv")
    df_IPR['RS'] = df_IPR['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})
    df_VPRG['RS'] = df_VPRG['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})
    df_IPR_TS3 = df_IPR[df_IPR['TS'] == 'TS3']
    df_IPR_TS4 = df_IPR[df_IPR['TS'] == 'TS4']
    plt.figure()
    sns.boxplot(x='RS',y='IPR', hue="wind", data=df_IPR_TS3, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.85)
    plt.legend(ncol=3)
    plt.figure()
    sns.boxplot(x='RS',y='IPR', hue="wind", data=df_IPR_TS4, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.85)
    plt.legend(ncol=3)
    plt.figure()
    sns.boxplot(x='RS',y='VPRG', hue="wind", data=df_VPRG, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.85)
    plt.legend(ncol=3)
    plt.show()

    return


#############################
# PLOTTERS FOR HYPOTHESIS 3 #
#############################

def determine_distance_boundaries_for_plot():
    # check windspeeds in conflict files
    distances = []
    for TS in [2, 4]:
        for RS in np.arange(1,5):
            conflict_report_json = open("thesis_tools/results/reports/conflict_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            conflict_data = json.load(conflict_report_json)
            conflict_report_json.close()
            for i in range(len(conflict_data)):
                distances.append(conflict_data[i]['dist_gf0'])
                distances.append(conflict_data[i]['dist_gf1'])
    distances.sort()
    first_boundary_idx = int(len(distances)/3) 
    second_boundary_idx = int(len(distances)//3 * 2)
    first_boundary = distances[first_boundary_idx]
    second_boundary_idx = distances[second_boundary_idx]
    max_distance = max(distances)
    return first_boundary, second_boundary_idx, max_distance

def IPR_VPRG_distance_box_whiskerplot_creator():
    df_IPR = pd.read_csv("thesis_tools/results/performance/hypothesis3_ipr.csv")
    df_VPRG = pd.read_csv("thesis_tools/results/performance/hypothesis3_vprg.csv")
    df_IPR['RS'] = df_IPR['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})
    df_VPRG['RS'] = df_VPRG['RS'].map({'RS1': 'OPT', 'RS2': 'DEST', 'RS3': 'HDG'})
    df_IPR_TS2 = df_IPR[df_IPR['TS'] == 'TS2']
    df_IPR_TS4 = df_IPR[df_IPR['TS'] == 'TS4']
    df_VPRG_TS2 = df_VPRG[df_VPRG['TS'] == 'TS2']
    df_VPRG_TS4 = df_VPRG[df_VPRG['TS'] == 'TS4']
    plt.figure()
    sns.boxplot(x='RS',y='IPR', hue="gf_dist", data=df_IPR_TS2,palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.85)
    plt.legend(ncol=3)
    plt.figure()
    sns.boxplot(x='RS',y='IPR', hue="gf_dist", data=df_IPR_TS4, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.85)
    plt.legend(ncol=3)
    plt.figure()
    sns.boxplot(x='RS',y='VPRG', hue="gf_dist", data=df_VPRG_TS2, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.85)
    plt.legend(ncol=3)
    plt.figure()
    sns.boxplot(x='RS',y='VPRG', hue="gf_dist", data=df_VPRG_TS4, palette="Greys")
    plt.grid(axis='y')
    plt.ylim(bottom=.85)
    plt.legend(ncol=3)
    plt.show()

    print(df_IPR)

    return

##############################
# Probability and statistics #
##############################

# hypothesis 1
def hypothesis1_data_generator():
    data_IPR = []
    for TS in np.arange(1,5):
        for RS in np.arange(1,4):
            IPRs = IPR_sample_calculator(TS, RS, sample_size)
            for k in range(len(IPRs)):
                data_IPR.append(['TS' + str(TS), 'RS' + str(RS), IPRs[k]])
    df_IPR = pd.DataFrame(data_IPR, columns = ["TS", "RS", "IPR"])
    df_IPR.to_csv("thesis_tools/results/performance/hypothesis1_ipr.csv")

    data_VPRG = []
    for TS in [2, 4]:
        for RS in np.arange(1,4):
            VPRGs = CVPRG_sample_calculator(TS, RS, sample_size)
            for k in range(len(VPRGs)):
                data_VPRG.append(['TS' + str(TS), 'RS' + str(RS), VPRGs[k]])
    df_VPRG = pd.DataFrame(data_VPRG, columns = ["TS", "RS", "VPRG"])
    df_VPRG.to_csv("thesis_tools/results/performance/hypothesis1_vprg.csv")

def Wilcoxon_hypothesis1():
    # Load data
    df_IPR = pd.read_csv("thesis_tools/results/performance/hypothesis1_ipr.csv")

    # check for each ruleset if TS1 vs TS2 and TS3 VS TS4
    p_data = []
    for TS in [[1,2],[3,4]]:
        for RS in np.arange(1,4):
            data = [[],[]]
            data[0] = df_IPR[df_IPR['TS'] == 'TS'  + str(TS[0])][df_IPR['RS'] == 'RS'  + str(RS)]["IPR"].tolist()
            data[1] = df_IPR[df_IPR['TS'] == 'TS'  + str(TS[1])][df_IPR['RS'] == 'RS'  + str(RS)]["IPR"].tolist()
            n = min(len(data[0]), len(data[1]))
            data[0] = random.sample(data[0], n)
            data[1] = random.sample(data[1], n)
            T, p = sp.stats.wilcoxon(data[0], data[1], zero_method="zsplit")
            p_data.append(['TS' + str(TS[0]), 'TS' + str(TS[1]), 'RS' + str(RS), 'IPR', n, p])
    df = pd.DataFrame(p_data, columns = ["TS_0", "TS_1", "RS", "variable", "n", "p-value"])
    df.to_csv("thesis_tools/results/pvalues/hypothesis1.csv")

# hypothesis 2
def hypothesis2_data_generator(first_slice, second_slice, sample_size):
    data_IPR = []
    data_VPRG = []
    for TS in [3, 4]:
        for RS in np.arange(1,4):
            grouped_data = [[], [], []]
            conflict_report_json = open("thesis_tools/results/reports/conflict_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            conflict_data = json.load(conflict_report_json)
            conflict_report_json.close()
            random.shuffle(conflict_data)
            for i in range(len(conflict_data)):
                # Sort conflict data
                if conflict_data[i]['windspeed'] < first_slice:
                    grouped_data[0].append(conflict_data[i])
                elif conflict_data[i]['windspeed'] < second_slice:
                    grouped_data[1].append(conflict_data[i])
                else:
                    grouped_data[2].append(conflict_data[i])
            
            # now sample data to obtain IPR
            for i in range(3):
                if i == 0:
                    wind_str = "low"
                elif i == 1:
                    wind_str = "medium"
                else:
                    wind_str = "strong"

                for j in range(int(len(grouped_data[i]) / sample_size)):
                    n_intrusion = 0
                    n_violations = 0
                    for k in range(sample_size):
                        idx = j * sample_size + k
                        n_intrusion += grouped_data[i][idx]['intrusion']
                        if TS == 4:
                            n_violations += grouped_data[i][idx]['violated_gf0']
                            n_violations += grouped_data[i][idx]['violated_gf1']
                            
                    if TS == 4:
                        VPRG = (2. * sample_size - n_violations) / (2. * sample_size)
                        data_VPRG.append(['TS' + str(TS), 'RS' + str(RS), wind_str, VPRG])
                    IPR = (sample_size - n_intrusion) / sample_size
                    data_IPR.append(['TS' + str(TS), 'RS' + str(RS), wind_str, IPR])

    df_IPR = pd.DataFrame(data_IPR, columns = ["TS", "RS", "wind", "IPR"])
    df_VPRG = pd.DataFrame(data_VPRG, columns = ["TS", "RS", "wind", "VPRG"])

    df_IPR.to_csv("thesis_tools/results/performance/hypothesis2_ipr.csv")
    df_VPRG.to_csv("thesis_tools/results/performance/hypothesis2_vprg.csv")

def Wilcoxon_hypothesis2():
    # Load data
    df_IPR = pd.read_csv("thesis_tools/results/performance/hypothesis2_ipr.csv")
    df_VPRG = pd.read_csv("thesis_tools/results/performance/hypothesis2_vprg.csv")

    # check for each windspeed/ruleset combination TS3 and TS4
    p_data = []
    for TS in [3, 4]:
        for RS in np.arange(1,4):
            for wind in [["low", "medium"], ["medium", "strong"]]:
                data = [[], []]
                data[0] = df_IPR[df_IPR['TS'] == 'TS' + str(TS)][df_IPR['RS'] == 'RS' + str(RS)][df_IPR['wind'] == wind[0]]["IPR"].tolist()
                data[1] = df_IPR[df_IPR['TS'] == 'TS' + str(TS)][df_IPR['RS'] == 'RS' + str(RS)][df_IPR['wind'] == wind[1]]["IPR"].tolist()
                n = min(len(data[0]), len(data[1]))
                data[0] = random.sample(data[0], n)
                data[1] = random.sample(data[1], n)
                T, p = sp.stats.wilcoxon(data[0], data[1], zero_method="zsplit")
                p_data.append(['TS' + str(TS), 'RS' + str(RS), wind[0], wind[1],'IPR', n, p])

    for TS in [4]:
        for RS in np.arange(1,4):
            for wind in [["low", "medium"], ["medium", "strong"]]:
                data = [[], []]
                data[0] = df_VPRG[df_VPRG['TS'] == 'TS' + str(TS)][df_VPRG['RS'] == 'RS' + str(RS)][df_VPRG['wind'] == wind[0]]["VPRG"].tolist()
                data[1] = df_VPRG[df_VPRG['TS'] == 'TS' + str(TS)][df_VPRG['RS'] == 'RS' + str(RS)][df_VPRG['wind'] == wind[1]]["VPRG"].tolist()
                n = min(len(data[0]), len(data[1]))
                data[0] = random.sample(data[0], n)
                data[1] = random.sample(data[1], n)
                T, p = sp.stats.wilcoxon(data[0], data[1], zero_method="zsplit")
                p_data.append(['TS' + str(TS), 'RS' + str(RS), wind[0], wind[1],'VPRG', n, p])

    df = pd.DataFrame(p_data, columns = ["TS", "RS", "wind_0", "wind_1", "variable", "n", "p-value"])
    df.to_csv("thesis_tools/results/pvalues/hypothesis2.csv")

# hypothesis 3
def hypothesis3_data_generator(first_slice, second_slice, sample_size):
    data_IPR = []
    data_VPRG = []
    for TS in [2, 4]:
        for RS in np.arange(1,4):
            grouped_data0 = [[], [], []]
            grouped_data1 = [[], [], []]
            conflict_report_json = open("thesis_tools/results/reports/conflict_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            conflict_data = json.load(conflict_report_json)
            conflict_report_json.close()
            random.shuffle(conflict_data)
            for i in range(len(conflict_data)):
                # Sort conflict data
                if conflict_data[i]['dist_gf0'] < first_slice:
                    grouped_data0[0].append(conflict_data[i])
                elif conflict_data[i]['dist_gf0'] < second_slice:
                    grouped_data0[1].append(conflict_data[i])
                else:
                    grouped_data0[2].append(conflict_data[i])

                if conflict_data[i]['dist_gf1'] < first_slice:
                    grouped_data1[0].append(conflict_data[i])
                elif conflict_data[i]['dist_gf1'] < second_slice:
                    grouped_data1[1].append(conflict_data[i])
                else:
                    grouped_data1[2].append(conflict_data[i])
            
            # now sample data to obtain IPR
            for i in range(3):
                if i == 0:
                    dist_str = "small"
                elif i == 1:
                    dist_str = "medium"
                else:
                    dist_str = "large"

                for j in range(int(len(grouped_data0[i]) / sample_size)):
                    n_intrusion = 0
                    n_violations = 0
                    for k in range(sample_size):
                        idx = j * sample_size + k
                        n_intrusion += grouped_data0[i][idx]['intrusion']
                        n_violations += grouped_data0[i][idx]['violated_gf0']
                    VPRG = (sample_size - n_violations) / sample_size
                    data_VPRG.append(['TS' + str(TS), 'RS' + str(RS), dist_str, VPRG])
                    IPR = (sample_size - n_intrusion) / sample_size
                    data_IPR.append(['TS' + str(TS), 'RS' + str(RS), dist_str, IPR])

                for j in range(int(len(grouped_data1[i]) / sample_size)):
                    n_intrusion = 0
                    n_violations = 0
                    for k in range(sample_size):
                        idx = j * sample_size + k
                        n_intrusion += grouped_data1[i][idx]['intrusion']
                        n_violations += grouped_data1[i][idx]['violated_gf1']
                    VPRG = (sample_size - n_violations) / sample_size
                    data_VPRG.append(['TS' + str(TS), 'RS' + str(RS), dist_str, VPRG])
                    IPR = (sample_size - n_intrusion) / sample_size
                    data_IPR.append(['TS' + str(TS), 'RS' + str(RS), dist_str, IPR])

    df_IPR = pd.DataFrame(data_IPR, columns = ["TS", "RS", "gf_dist", "IPR"])
    df_VPRG = pd.DataFrame(data_VPRG, columns = ["TS", "RS", "gf_dist", "VPRG"])
    df_IPR.to_csv("thesis_tools/results/performance/hypothesis3_ipr.csv")
    df_VPRG.to_csv("thesis_tools/results/performance/hypothesis3_vprg.csv")

def Wilcoxon_hypothesis3():
    # Load data
    df_IPR = pd.read_csv("thesis_tools/results/performance/hypothesis3_ipr.csv")
    df_VPRG = pd.read_csv("thesis_tools/results/performance/hypothesis3_vprg.csv")

    # check for each geofence distance/ruleset combination TS3 and TS4
    p_data = []
    for TS in [2, 4]:
        for RS in np.arange(1,4):
            for dist in [["small", "medium"], ["medium", "large"]]:
                data = [[], []]
                data[0] = df_IPR[df_IPR['TS'] == 'TS' + str(TS)][df_IPR['RS'] == 'RS' + str(RS)][df_IPR['gf_dist'] == dist[0]]["IPR"].tolist()
                data[1] = df_IPR[df_IPR['TS'] == 'TS' + str(TS)][df_IPR['RS'] == 'RS' + str(RS)][df_IPR['gf_dist'] == dist[1]]["IPR"].tolist()
                n = min(len(data[0]), len(data[1]))
                data[0] = random.sample(data[0], n)
                data[1] = random.sample(data[1], n)
                T, p = sp.stats.wilcoxon(data[0], data[1], zero_method="zsplit")
                p_data.append(['TS' + str(TS), 'RS' + str(RS), dist[0], dist[1],'IPR', n, p])

    for TS in [4]:
        for RS in np.arange(1,4):
            for dist in [["small", "medium"], ["medium", "large"]]:
                data = [[], []]
                data[0] = df_VPRG[df_VPRG['TS'] == 'TS' + str(TS)][df_VPRG['RS'] == 'RS' + str(RS)][df_VPRG['gf_dist'] == dist[0]]["VPRG"].tolist()
                data[1] = df_VPRG[df_VPRG['TS'] == 'TS' + str(TS)][df_VPRG['RS'] == 'RS' + str(RS)][df_VPRG['gf_dist'] == dist[1]]["VPRG"].tolist()
                n = min(len(data[0]), len(data[1]))
                data[0] = random.sample(data[0], n)
                data[1] = random.sample(data[1], n)
                T, p = sp.stats.wilcoxon(data[0], data[1], zero_method="zsplit")
                p_data.append(['TS' + str(TS), 'RS' + str(RS), dist[0], dist[1],'VPRG', n, p])

    df = pd.DataFrame(p_data, columns = ["TS", "RS", "dist_0", "dist_1", "variable", "n", "p-value"])
    df.to_csv("thesis_tools/results/pvalues/hypothesis3.csv")

#hypothesis1_data_generator()
#hypothesis2_data_generator(5, 10, sample_size)
#hypothesis3_data_generator(150, 250, sample_size)

#Wilcoxon_hypothesis1()
#Wilcoxon_hypothesis2()
#Wilcoxon_hypothesis3()

#IPR_box_whiskerplot_creator()
#CVPRG_box_whiskerplot_creator()
#print(determine_wind_boundaries_for_plot())
#IPR_VPRG_wind_box_whiskerplot_creator()
#print(determine_distance_boundaries_for_plot())
IPR_VPRG_distance_box_whiskerplot_creator()