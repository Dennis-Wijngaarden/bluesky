import numpy as np 
import matplotlib.pyplot as plt 
import os
import json
import multiprocessing

import warnings
warnings.simplefilter("ignore", UserWarning)

log_interval = 1.0 # [s]
R_pz = 50.

def analyse_reference(TS, path):
    # Read files in path folder
    fllog_files = np.array(os.listdir(path + "/FLLOGS"))
    gflog_files = np.array(os.listdir(path + "/GFLOGS"))

    indices = [i for i, v in enumerate(fllog_files) if 'TS' + str(TS) in v]
    
    # Sort filenames by log types in the following lists
    fllog_files = fllog_files[indices]
    gflog_files = gflog_files[indices]
    
    # Store data in arrays using genfromtxt
    fllog_raw_data = []
    fllog_callsigns = []
    gflog_raw_data = []
    gflog_callsigns = []

    invalid_missions = []
    csv_text = "#time0 [s], distance0 [-], leg_finished_0 [-], min_gf_dist0 [m], violation0 [-], time1 [s], distance1 [-], leg_finished_1 [-], min_gf_dist1 [m], violation1 [-], scen_valid [-]\n"

    for i in range(len(fllog_files)):
        print("analysing reference TS"+ str(TS) + " number: " + str(i))
        fllog_raw_data = np.genfromtxt(path + '/FLLOGS/' + fllog_files[i], delimiter = ',', usecols = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        fllog_callsigns = np.genfromtxt(path + '/FLLOGS/' + fllog_files[i], delimiter = ',', usecols = [1], dtype=str)

        gflog_raw_data = np.genfromtxt(path + '/GFLOGS/' + gflog_files[i], delimiter = ',', usecols = [3])
        gflog_callsigns = np.genfromtxt(path + '/GFLOGS/' + gflog_files[i], delimiter = ',', usecols = [1], dtype=str)

        # First leg data
        # check if first leg is finished
        time0, distance0, leg_finished0 = filter_flight_reference_data(fllog_raw_data, fllog_callsigns, "UAV0")
        time1, distance1, leg_finished1 = filter_flight_reference_data(fllog_raw_data, fllog_callsigns, "UAV1")

        # Geofence data
        # check if negative geofence distance
        min_gf_dist0, violation0 = filter_geofence_reference_date(gflog_raw_data, gflog_callsigns, "UAV0")
        min_gf_dist1, violation1 = filter_geofence_reference_date(gflog_raw_data, gflog_callsigns, "UAV1")

        # validaty data
        validity = leg_finished0 and leg_finished1 and not violation0 and not violation1

        # write line to csv
        csv_text += str(time0) + ', ' + str(distance0) + ', ' + str(leg_finished0) + ', ' + str(min_gf_dist0) + ", " + str(violation0) + ', ' +\
                    str(time1) + ', ' + str(distance1) + ', ' + str(leg_finished1) + ', ' + str(min_gf_dist1) + ", " + str(violation1) + ', ' +\
                    str(validity) + '\n'

        if not validity:
            invalid_missions.append(i)
    
    # write reference csv to file
    reference_csv_file = open('thesis_tools/results/reference/ref_TS' + str(TS) + '.csv', 'w')
    reference_csv_file.write(csv_text)
    reference_csv_file.close()

    return invalid_missions


def filter_flight_reference_data(raw_fl_data_log, callsigns_log, callsign):
    indices = np.where(callsigns_log == callsign)
    raw_fl_data = raw_fl_data_log[indices]

    # check line where active waypoint changes
    lat0 = raw_fl_data[0][8]

    try:
        index = np.where(raw_fl_data[:,8] != lat0)[0][0]
        leg_finished = True
        time = raw_fl_data[index][0]
        distance = raw_fl_data[index][3]
    except:
        leg_finished = False
        time = 0.
        distance = 0.
    
    return time, distance, leg_finished

def filter_geofence_reference_date(raw_gf_data_log, callsigns_log, callsign):
    # Determine rows for which 'callsign' is the ownship
    indices = np.where(callsigns_log == callsign)
    raw_gf_data = raw_gf_data_log[indices]
    min_gf_dist = min(raw_gf_data)
    violation = min_gf_dist < 0
    return min_gf_dist, violation

def analyse_ruleset_in_testseries(TS, RS, geofence_defined, wind_defined, invalid_indices, path):
    # Read files in path folder
    fllog_files = np.array(os.listdir(path + "/FLLOGS"))
    conflog_files = np.array(os.listdir(path + "/CONFLOGS"))
    gflog_files = np.array(os.listdir(path + "/GFLOGS"))

    file_indices = [i for i, v in enumerate(fllog_files) if ((('TS' + str(TS)) in v) and (('RS' + str(RS)) in v))]
    
    # Sort filenames by log types in the following lists
    fllog_files = fllog_files[file_indices]
    conflog_files = conflog_files[file_indices]
    gflog_files = gflog_files[file_indices]

    # Delete the invalid_indices
    valid_indices = list(set(range(len(file_indices))) - set(invalid_indices))
    fllog_files = fllog_files[valid_indices]
    conflog_files = conflog_files[valid_indices]
    gflog_files = gflog_files[valid_indices]

    # Load wind data id defined
    if wind_defined:
        # open wind json file
        # load route data from route.json
        wind_json = open("thesis_tools/data/wind.json", "r")
        wind_data = json.load(wind_json)
        wind_json.close()

    # Define reports
    conflict_report = []
    scenario_report = []

    for i in range(len(valid_indices)):
        # status print
        print("Analysing TS" + str(TS) + " RS" + str(RS) + " number: " + str(valid_indices[i]))
        fllog_raw_data = np.genfromtxt(path + '/FLLOGS/' + fllog_files[i], delimiter = ',', usecols = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        fllog_callsigns = np.genfromtxt(path + '/FLLOGS/' + fllog_files[i], delimiter = ',', usecols = [1], dtype=str)

        conflog_raw_data = np.genfromtxt(path + '/CONFLOGS/' + conflog_files[i], delimiter = ',', usecols = [0, 3, 4, 5, 6, 7])
        conflog_callsigns = np.genfromtxt(path + '/CONFLOGS/' + conflog_files[i], delimiter = ',', usecols = [1, 2], dtype=str)
    	
        if geofence_defined:
            gflog_raw_data = np.genfromtxt(path + '/GFLOGS/' + gflog_files[i], delimiter = ',', usecols = [0, 3])
            gflog_callsigns = np.genfromtxt(path + '/GFLOGS/' + gflog_files[i], delimiter = ',', usecols = [1], dtype=str)

        if wind_defined:
            windspeed = wind_data[valid_indices[i]]['speed']

        # Analyse conflicts
        t_start_conflict, t_conflict, min_dist_conflict, min_t_los_conflict = get_conflict_variables(conflog_raw_data, conflog_callsigns, 'UAV0')

        # Analyse geofences
        if geofence_defined:
            # Get some extra variables regarding geofences
            dist_geo_0, geo_violated0 = get_geo_variables_for_conflict(t_start_conflict, gflog_raw_data, gflog_callsigns, 'UAV0')
            dist_geo_1, geo_violated1 = get_geo_variables_for_conflict(t_start_conflict, gflog_raw_data, gflog_callsigns, 'UAV1')

            # min distance encountered wrt the geofence
            min_dist_geofence0, min_dist_geofence_time0 = get_geofence_variables(gflog_raw_data, gflog_callsigns, 'UAV0')
            min_dist_geofence1, min_dist_geofence_time1 = get_geofence_variables(gflog_raw_data, gflog_callsigns, 'UAV1')

        # write conflict report
        for j in range(len(t_start_conflict)):
            conflict_line = {} # dictionary with conflict variable
            conflict_line['t_conflict'] = t_conflict[j]
            conflict_line['intrusion'] = int(min_dist_conflict[j] < R_pz)
            conflict_line['scenario'] = valid_indices[i]
            if wind_defined:
                conflict_line['windspeed'] = windspeed
            if geofence_defined:
                conflict_line['dist_gf0'] = dist_geo_0[j]
                conflict_line['dist_gf1'] = dist_geo_1[j]
                conflict_line['violated_gf0'] = int(geo_violated0[j])
                conflict_line['violated_gf1'] = int(geo_violated1[j])
            conflict_report.append(conflict_line)
        
        # write scenario report
        scenario_line = {} # dictionary with scenario variables
        scenario_line['scenario'] = valid_indices[i]
        scenario_line["n_conflicts"] = len(t_start_conflict)
        scenario_line['n_PZ_violated'] = 0
        if len(min_dist_conflict) > 0:
            scenario_line['min_rel_dist'] = min(min_dist_conflict)
            for j in range(len(min_dist_conflict)):
                if (min_dist_conflict[j] < R_pz):
                    scenario_line['n_PZ_violated'] += 1

        if geofence_defined:
            scenario_line['min_dist_gf0'] = min_dist_geofence0
            scenario_line['min_dist_gf1'] = min_dist_geofence1
            scenario_line['gf_violated0'] = int(min_dist_geofence0 < 0)
            scenario_line['gf_violated1'] = int(min_dist_geofence1 < 0)

        if wind_defined:
            scenario_line['windspeed'] = windspeed
        scenario_report.append(scenario_line)
        

    location = "thesis_tools/results/reports"
    conflict_report_json = json.dumps(conflict_report, indent=4)
    conflict_report_file = open(location + "/conflict_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "w")
    conflict_report_file.write(conflict_report_json)
    conflict_report_file.close()

    location = "thesis_tools/results/reports"
    scenario_report_json = json.dumps(scenario_report, indent=4)
    scenario_report_file = open(location + "/scenario_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "w")
    scenario_report_file.write(scenario_report_json)
    scenario_report_file.close()

    return

def filter_raw_conflict_data(raw_conf_data_log, callsigns_log, callsign):
    # Determine rows for which 'callsign' is the ownship
    if (len(callsigns_log) > 0):
        indices = np.where(callsigns_log[:,0] == callsign)
    else:
        indices = np.array([], dtype = int)
    raw_conflict_data = raw_conf_data_log[indices]

    # Now determine start and end index of each conflict
    conflicts_detected = False
    if (len(raw_conflict_data) > 0):
        conflicts_detected = True
        dt = np.diff(raw_conflict_data[:,0])
    
    conflict_data = []
    if (conflicts_detected):
        if (len(dt) == 0):
            conflict_data.append(raw_conflict_data)
        else:
            for i in range(len(dt)):
                conflict_indices = np.where(dt > log_interval)[0] + 1
                conflict_indices = np.append(np.array([0]), conflict_indices)
                conflict_indices = np.append(conflict_indices, len(raw_conflict_data))
            for i in range(len(conflict_indices) - 1):
                conflict_data.append(raw_conflict_data[conflict_indices[i] : conflict_indices[i + 1]])
    
    return conflict_data

def filter_raw_geofence_data(raw_gf_data_log, callsigns_log, callsign):
    # Determine rows for which 'callsign' is the ownship
    indices = np.where(callsigns_log == callsign)
    geofence_data = raw_gf_data_log[indices]
    
    return geofence_data

def get_conflict_variables(conflog_raw_data, conflog_callsigns, callsign):
    # Determine each single conflict for each vehicle
    t_start_conflict = [] # Start time of the conflict
    t_conflict = [] # Duration of the conflict
    min_dist_conflict = [] # Minimal dist encountered during conflict
    min_t_los_conflict = [] # Minimal time to los during conflict

    conflict_data = filter_raw_conflict_data(conflog_raw_data, conflog_callsigns, callsign)

    for i in range(len(conflict_data)):
        t_start_conflict.append(conflict_data[i][0][0])
        t_conflict.append(len(conflict_data[i]) * log_interval)

        dist = []
        t_los = []
        for j in range(len(conflict_data[i])):
            dist.append(conflict_data[i][j][2])
            t_los.append(conflict_data[i][j][5])
        min_dist_conflict.append(min(dist))
        min_t_los_conflict.append(min(t_los))
    return t_start_conflict, t_conflict, min_dist_conflict, min_t_los_conflict

def get_geo_variables_for_conflict(t_start_conflict, gflog_raw_data, gflog_callsigns, callsign):
    # Create lists with distances to geofence and geofence violated after conflict
    dist_geo = []
    geo_violated = []

    # Determine for which rows correspond to the ownship
    gf_data = np.array(filter_raw_geofence_data(gflog_raw_data, gflog_callsigns, callsign))

    # Check if there is a violation in general
    gf_violation_general = min(gf_data[:,1]) < 0
    for i in range(len(t_start_conflict)):
        idx_start_conflict = np.where(gf_data[:,0] == t_start_conflict[i])[0][0]
        # Find distance w.r.t. geofence at start of conflict
        dist_geo.append(gf_data[idx_start_conflict][1])
        # Check if geofence is violated after conflict (after current conflict, but before next conflict or end)
        if (gf_violation_general):
            # First check if vehicle is at begin of leg inside geofence
            in_geofence = gf_data[idx_start_conflict][1] > 0
            if in_geofence:
                if (i == (len(t_start_conflict) - 1)):
                    # This conflict is the final conflict
                    gf_violated = min(gf_data[idx_start_conflict:,1]) < 0
                else:
                    # This conflict is not the final conflict
                    idx_start_conflict_next = np.where(gf_data[:,0] == t_start_conflict[i+1])[0][0]
                    gf_violated = min(gf_data[idx_start_conflict:idx_start_conflict_next,1]) < 0
            else:
                # Geofence not violated due to conflict when already outside geofence at begin of conflict
                gf_violated = False
            geo_violated.append(gf_violated)
        else:
            geo_violated.append(False)
            
    return dist_geo, geo_violated

def get_geofence_variables(gflog_raw_data, gflog_callsigns, callsign):
    # Determine for which rows correspond to the ownship
    gf_data = filter_raw_geofence_data(gflog_raw_data, gflog_callsigns, callsign)
    rel_distances = np.array(gf_data)[:,1].tolist()
    min_dist_geofence = min(rel_distances)

    min_dist_geofence_time = np.where(gf_data[:,1] == min_dist_geofence)[0][0]
    
    return min_dist_geofence, min_dist_geofence_time

def get_geofence_conflict_variables(conflog_raw_data, conflog_callsigns, gflog_raw_data, gflog_callsigns, callsign):
    # distance wrt to geodence at start of conflict
    dist_geo_start = [] # distance wrt geofence at the start of a conflict
    dist_geo_min = [] # minimum distance wrt gefoence during conflict
    time_geo_min = [] # Time at which distacne wrt geofence is at minimum

    for i in range(len(conflog_raw_data)):
        # filter conflict data
        conflict_data = filter_raw_conflict_data(conflog_raw_data[i], conflog_callsigns[i], callsign)
        gf_data = filter_raw_geofence_data(gflog_raw_data[i], gflog_callsigns[i], callsign)
        indices_geo_start = [] # indices of gf_data where conflicts start
        dist_geo_start.append([])
        dist_geo_min.append([])
        time_geo_min.append([])
        for j in range(len(conflict_data)):
            indices_geo_start.append(np.where(gf_data[:,0] == conflict_data[j][0][0])[0][0])
            dist_geo_start[i].append(gf_data[indices_geo_start[j]][1])
            distances_geo = []
            for k in range(len(conflict_data[j])):
                index_gf_data = np.where(gf_data[:,0] == conflict_data[j][k][0])
                distances_geo.append(gf_data[index_gf_data[0][0]][1])
            dist_geo_min[i].append(min(distances_geo))
    return dist_geo_start, dist_geo_min

def report_IRPZ(TS_list, RS_list):
    for TS in TS_list:
        for RS in RS_list:
            # open json file:
            simple_report_json = open("thesis_tools/results/reports/simple_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            report_data = json.load(simple_report_json)
            simple_report_json.close()

            # Now determine irpz
            n_sim = len(report_data)
            n_los = 0
            for i in range(n_sim):
                if report_data[i]['n_PZ_violated'] > 0:
                    n_los += 1
            print("IRPZ TS" + str(TS) + " RS" + str(RS) + ": " + str(n_los / n_sim))
    return

def report_IPR(TS_list, RS_list):
    for TS in TS_list:
        for RS in RS_list:
            # open json file:
            simple_report_json = open("thesis_tools/results/reports/simple_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            report_data = json.load(simple_report_json)
            simple_report_json.close()

            # Now determine IPR
            n_sim = len(report_data)
            n_cfl = 0
            n_los = 0
            for i in range(n_sim):
                n_cfl += report_data[i]['n_conflicts']
                n_los += report_data[i]['n_PZ_violated']
            print("IPR TS" + str(TS) + " RS" + str(RS) + ": " + str((n_cfl - n_los) / n_cfl))
    return

def report_VRG(TS_list, RS_list):
    for TS in TS_list:
        for RS in RS_list:
            # open json file:
            simple_report_json = open("thesis_tools/results/reports/simple_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "r")
            report_data = json.load(simple_report_json)
            simple_report_json.close()

            # Now determine irpz
            n_sim = len(report_data)
            n_gv = 0
            for i in range(n_sim):
                if report_data[i]["gf_violated0"] == 1:
                    n_gv += 1
                if report_data[i]["gf_violated1"] == 1:
                    n_gv += 1
            print("VRG TS" + str(TS) + " RS" + str(RS) + ": " + str(n_gv / (n_sim * 2.)))

    return

if __name__ == '__main__':
    invalid_indices0 = analyse_reference(0, "output")
    invalid_indices5 = analyse_reference(5, "output")
    invalid_indices = list(set(invalid_indices0) | set(invalid_indices5)) 
    invalid_indices.sort()

    TS_list = [1,2,3,4]
    RS_list = [1,2,3]
    geo_list = [False, True, False, True]
    wind_list = [False, False, True, True]

    idx = 0
    processes = []
    for TS in TS_list:
        for RS in RS_list:
            processes.append(multiprocessing.Process(target = analyse_ruleset_in_testseries, args = (TS, RS, geo_list[idx], wind_list[idx], invalid_indices, "output")))
            processes[-1].start()
        idx += 1