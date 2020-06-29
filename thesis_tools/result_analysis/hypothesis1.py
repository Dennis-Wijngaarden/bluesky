import numpy as np
import json

def IPR_analysis():
    # first come up with statistics and print
    
    for TS in [1,2,3,4]:
        for RS in [1,2,3]:
            conflict_report_json = open("thesis_tools/results/reports/conflict_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            conflict_data = json.load(conflict_report_json)
            conflict_report_json.close()

            scenario_report_json = open("thesis_tools/results/reports/scenario_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            scenario_data = json.load(scenario_report_json)
            scenario_report_json.close()

            # print statistics
            print('\nReport TS' + str(TS) + ' RS' + str(RS) + ' :')
            n_conflicts = len(conflict_data)
            n_LoS = 0
            for i in range(len(conflict_data)):
                n_LoS += conflict_data[i]['intrusion']
            scenario_list = []
            n_conflicts_scenario_list = []
            for i in range(len(scenario_data)):
                scenario_list.append(scenario_data[i]['scenario'])
                n_conflicts_scenario_list.append(scenario_data[i]['n_conflicts']) 
            n_scenarios = len(scenario_list)
            n_conflicts_avg = n_conflicts / n_scenarios
            n_conflicts_med = np.median(n_conflicts_scenario_list)
            print('n_conflicts: ' + str(n_conflicts))
            print('n_LoS: ' + str(n_LoS))
            print('IPR: ' + str((n_conflicts - n_LoS) / n_conflicts))
            print('n_scenarios: ' + str(n_scenarios))
            print('n_conflicts_avg: ' + str(n_conflicts_avg))
            print('n_conflicts_med: ' + str(int(n_conflicts_med)))

# Find out why DEST performs worst in terms of IPR for all test series
def DEST_IPR():
    # Find indices that succeed in OPT and HDG but not in DEST
    test_numbers = []
    for TS in [1,2,3,4]:
        test_numbers.append([])
        scenario_report_DEST_json = open("thesis_tools/results/reports/scenario_report_TS" + str(TS) + "_RS" + str(2) + ".json", "r")
        scenario_data_DEST = json.load(scenario_report_DEST_json)
        scenario_report_DEST_json.close()

        scenario_data = []
        for RS in [1,3]:
            scenario_data.append([])
            scenario_report_json = open("thesis_tools/results/reports/scenario_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            scenario_data[-1] = json.load(scenario_report_json)
            scenario_report_json.close()


        for i in range(len(scenario_data_DEST)):
            violation = bool(scenario_data[0][i]['n_PZ_violated'] or scenario_data[1][i]['n_PZ_violated'])
            if (scenario_data_DEST[i]['n_PZ_violated'] and not violation):
                test_numbers[-1].append(int(scenario_data_DEST[i]['scenario']))


    return test_numbers

#IPR_analysis()
print(DEST_IPR()[0])
