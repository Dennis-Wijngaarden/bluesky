import numpy as np 
import matplotlib.pyplot as plt 
import os
import json

log_interval = 1.0 # [s]
R_pz = 25.

def analyse_reference(TS, path):
    # Read files in path folder
    files = os.listdir(path)
    
    # Sort filenames by log types in the following lists
    fllog_files = []
    gflog_files = []
    for i in range(len(files)):
        if (("TS" + str(TS)) in files[i]):
            if ("FLLOG" in files[i]):
                fllog_files.append(files[i])
            elif ("GFLOG" in files[i]):
                gflog_files.append(files[i])
    
    # Store data in arrays using genfromtxt
    fllog_raw_data = []
    fllog_callsigns = []
    gflog_raw_data = []
    gflog_callsigns = []

    for i in range(len(fllog_files)):
        fllog_raw_data.append(np.genfromtxt(path + '/' + fllog_files[i], delimiter = ',', usecols = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10]))
        fllog_callsigns.append(np.genfromtxt(path + '/' + fllog_files[i], delimiter = ',', usecols = [1], dtype=str))

        gflog_raw_data.append(np.genfromtxt(path + '/' + gflog_files[i], delimiter = ',', usecols = [3]))
        gflog_callsigns.append(np.genfromtxt(path + '/' + gflog_files[i], delimiter = ',', usecols = [1], dtype=str))

    invalid_missions = []
    csv_text = "#time0 [s], distance0 [-], leg_finished_0 [-], min_gf_dist0 [m], violation0 [-], time1 [s], distance1 [-], leg_finished_1 [-], min_gf_dist1 [m], violation1 [-], scen_valid [-]\n"
    for i in range(len(fllog_raw_data)):
        # First leg data
        # check if first leg is finished
        time0, distance0, leg_finished0 = filter_flight_reference_data(fllog_raw_data[i], fllog_callsigns[i], "UAV0")
        time1, distance1, leg_finished1 = filter_flight_reference_data(fllog_raw_data[i], fllog_callsigns[i], "UAV1")

        # Geofence data
        # check if negative geofence distance
        min_gf_dist0, violation0 = filter_geofence_reference_date(gflog_raw_data[i], gflog_callsigns[i], "UAV0")
        min_gf_dist1, violation1 = filter_geofence_reference_date(gflog_raw_data[i], gflog_callsigns[i], "UAV1")

        # validaty data
        validity = leg_finished0 and leg_finished1 and not violation0 and not violation1

        if not validity:
            invalid_missions.append(i)

        # write line to csv
        csv_text += str(time0) + ', ' + str(distance0) + ', ' + str(leg_finished0) + ', ' + str(min_gf_dist0) + ", " + str(violation0) + ', ' +\
                    str(time1) + ', ' + str(distance1) + ', ' + str(leg_finished1) + ', ' + str(min_gf_dist1) + ", " + str(violation1) + ', ' +\
                    str(validity) + '\n'
    
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
    lon0 = raw_fl_data[0][9]
    leg_finished = False
    time = 0.
    distance = 0.
    for i in range(len(raw_fl_data)):
        if ((raw_fl_data[i][8] != lat0) or (raw_fl_data[i][9] != lon0)):
            time = raw_fl_data[i][0]
            distance = raw_fl_data[i][3]
            leg_finished = True
            break
    return time, distance, leg_finished

def filter_geofence_reference_date(raw_gf_data_log, callsigns_log, callsign):
    # Determine rows for which 'callsign' is the ownship
    indices = np.where(callsigns_log == callsign)
    raw_gf_data = raw_gf_data_log[indices]
    min_gf_dist = min(raw_gf_data)
    violation = min_gf_dist < 0
    return min_gf_dist, violation

def analyse_ruleset_in_testseries(TS, RS, geofence_defined, path):
    # Read files in path folder
    files = os.listdir(path)

    # Sort filenames by log types in the following lists
    fllog_files = []
    conflog_files = []
    gflog_files = []
    for i in range(len(files)):
        if ((("TS" + str(TS)) in files[i]) and (("RS" + str(RS)) in files[i])):
            if ("FLLOG" in files[i]):
                fllog_files.append(files[i])
            elif ("CONFLOG" in files[i]):
                conflog_files.append(files[i])
            elif ("GFLOG" in files[i]):
                gflog_files.append(files[i])
    
    # Store data in arrays using genfromtxt
    fllog_raw_data = []
    fllog_callsigns = []
    conflog_raw_data = []
    conflog_callsigns = []
    gflog_raw_data = []
    gflog_callsigns = []

    for i in range(len(conflog_files)):
        fllog_raw_data.append(np.genfromtxt(path + '/' + fllog_files[i], delimiter = ',', usecols = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10]))
        fllog_callsigns.append(np.genfromtxt(path + '/' + fllog_files[i], delimiter = ',', usecols = [1], dtype=str))

        conflog_raw_data.append(np.genfromtxt(path + '/' + conflog_files[i], delimiter = ',', usecols = [0, 3, 4, 5, 6, 7]))
        conflog_callsigns.append(np.genfromtxt(path + '/' + conflog_files[i], delimiter = ',', usecols = [1, 2], dtype=str))
    	
        if geofence_defined:
            gflog_raw_data.append(np.genfromtxt(path + '/' + gflog_files[i], delimiter = ',', usecols = [0, 3]))
            gflog_callsigns.append(np.genfromtxt(path + '/' + gflog_files[i], delimiter = ',', usecols = [1], dtype=str))

    # Analyse conflicts
    t_start_conflict, t_conflict, min_dist_conflict, min_t_los_conflict = get_conflict_variables(conflog_raw_data, conflog_callsigns, 'UAV0')

    # Analyse geofences
    if geofence_defined:
        # min distance encountered wrt the geofence
        min_dist_geofence0, min_dist_geofence_time0 = get_geofence_variables(gflog_raw_data, gflog_callsigns, 'UAV0')
        min_dist_geofence1, min_dist_geofence_time1 = get_geofence_variables(gflog_raw_data, gflog_callsigns, 'UAV1')
        
        # geo/conflict variables
        dist_geo_start0, dist_geo_min0 = get_geofence_conflict_variables(conflog_raw_data, conflog_callsigns, gflog_raw_data, gflog_callsigns, 'UAV0')
        dist_geo_start1, dist_geo_min1 = get_geofence_conflict_variables(conflog_raw_data, conflog_callsigns, gflog_raw_data, gflog_callsigns, 'UAV1')

    # Determine set of valid scenarios
    invalid_indices0 = analyse_reference(0, "output") # Wind calm
    invalid_indices5 = analyse_reference(5, "output") # windy

    valid_indices = list(set(range(len(fllog_files))) - set(invalid_indices0) - set(invalid_indices5))

    # Now write imple json report
    simple_report = []
    for i in valid_indices:
        scenario_line = {} # dictionary with scenario variables
        scenario_line['scenario'] = i
        scenario_line["n_conflicts"] = len(min_dist_conflict[i])
        scenario_line['n_PZ_violated'] = 0
        if len(min_dist_conflict[i]) > 0:
            scenario_line['min_rel_dist'] = min(min_dist_conflict[i])
            for j in range(len(min_dist_conflict[i])):
                if (min_dist_conflict[i][j] < R_pz):
                    scenario_line['n_PZ_violated'] += 1

        if geofence_defined:
            scenario_line['min_dist_gf0'] = min_dist_geofence0[i]
            scenario_line['min_dist_gf1'] = min_dist_geofence1[i]
            if len(min_dist_conflict[i]) > 0:
                scenario_line['dist_geo_start0'] = min(dist_geo_start0[i])
                scenario_line['dist_geo_start1'] = min(dist_geo_start1[i])
            scenario_line['gf_violated0'] = int(min_dist_geofence0[i] < 0)
            scenario_line['gf_violated1'] = int(min_dist_geofence1[i] < 0)
        simple_report.append(scenario_line)

    location = "thesis_tools/results/reports"
    simple_report_json = json.dumps(simple_report, indent=4)
    simple_report_file = open(location + "/simple_report_TS" + str(TS) + "_RS" + str(RS) + ".json", "w")
    simple_report_file.write(simple_report_json)
    simple_report_file.close()

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
    for i in range(len(conflog_raw_data)):
        # Determine rows for which UAv0 is the ownship
        conflict_data = filter_raw_conflict_data(conflog_raw_data[i], conflog_callsigns[i], callsign)
        # Analyse each single conflict occured in scenario
        t_start_conflict.append([])
        t_conflict.append([])
        min_dist_conflict.append([])
        min_t_los_conflict.append([])
        for j in range(len(conflict_data)):
            t_start_conflict[i].append(conflict_data[j][0][0])
            t_conflict[i].append(len(conflict_data[j]) * log_interval)

            dist = []
            t_los = []
            for k in range(len(conflict_data[j])):
                dist.append(conflict_data[j][k][2])
                t_los.append(conflict_data[j][k][5])
            min_dist_conflict[i].append(min(dist))
            min_t_los_conflict[i].append(min(t_los))
    return t_start_conflict, t_conflict, min_dist_conflict, min_t_los_conflict

def get_geofence_variables(gflog_raw_data, gflog_callsigns, callsign):
    # Determube minimum geofence distance and violation
    min_dist_geofence = []
    min_dist_geofence_time = []
    for i in range(len(gflog_raw_data)):
        # Determine for which rows correspond to the ownship
        gf_data = filter_raw_geofence_data(gflog_raw_data[i], gflog_callsigns[i], callsign)
        rel_distances = np.array(gf_data)[:,1].tolist()
        min_dist_geofence.append(min(rel_distances))

        min_dist_geofence_time.append(np.where(gf_data[:,1] == min_dist_geofence[-1])[0][0])
    
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

def generate_reports(TS_list, RS_list, geo_list):
    geo_idx = 0
    for TS in TS_list:
        for RS in RS_list:
            print("Generate report TS" + str(TS) + " RS" + str(RS))
            analyse_ruleset_in_testseries(TS, RS, geo_list[geo_idx], "output")
        geo_idx += 1
    return

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

generate_reports([1,2,3,4], [1,2,3,4], [False, True, False, True])
report_IRPZ([1,2,3,4], [1,2,3,4])
report_IPR([1,2,3,4], [1,2,3,4])
report_VRG([2,4], [1,2,3,4])
#generate_reports([1,2,3,4], [1,2,3,4], [False, True, False, True])
#analyse_ruleset_in_testseries(2, 3, True, 'output')

#generate_performance_reports([0], [0])
#t1, dist1, t_los1, dist_geo1_0, dist_geo1_1 = analyse_ruleset_in_testseries(1, 1, "output")

#t1, dist1, t_los1, dist_geo1_0, dist_geo1_1 = analyse_ruleset_in_testseries(1, 1, "output")
#t3, dist3, t_los3, dist_geo3_0, dist_geo3_1 = analyse_ruleset_in_testseries(1, 3, "output")

#flat_dist1 = []
#flat_indices1 = []
#for i in range(len(dist1)):
#    for j in range(len(dist1[i])):
#        flat_dist1.append(dist1[i][j])
#        flat_indices1.append(i)
#
#flat_dist3 = []
#flat_indices3 = []
#for i in range(len(dist3)):
#    for j in range(len(dist3[i])):
#        flat_dist3.append(dist3[i][j])
#        flat_indices3.append(i)
#
#print(flat_indices1)
#plt.scatter(flat_indices3, flat_dist3)
#plt.scatter(flat_indices1, flat_dist1)
#plt.show()