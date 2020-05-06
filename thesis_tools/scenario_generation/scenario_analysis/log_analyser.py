import numpy as np 
import matplotlib.pyplot as plt 
import os

log_interval = 1.0 # [s]

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

        # write line to csv
        csv_text += str(time0) + ', ' + str(distance0) + ', ' + str(leg_finished0) + ', ' + str(min_gf_dist0) + ", " + str(violation0) + ', ' +\
                    str(time1) + ', ' + str(distance1) + ', ' + str(leg_finished1) + ', ' + str(min_gf_dist1) + ", " + str(violation1) + ', ' +\
                    str(validity) + '\n'
    
    # write reference csv to file
    reference_csv_file = open('thesis_tools/results/reference/ref_TS' + str(TS) + '.csv', 'w')
    reference_csv_file.write(csv_text)
    reference_csv_file.close()


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

analyse_reference(5, "output")

def analyse_ruleset_in_testseries(TS, RS, path):
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
    #gflog_raw_data = []
    #gflog_callsigns = []
    for i in range(len(conflog_files)):
        fllog_raw_data.append(np.genfromtxt(path + '/' + fllog_files[i], delimiter = ',', usecols = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10]))
        fllog_callsigns.append(np.genfromtxt(path + '/' + fllog_files[i], delimiter = ',', usecols = [1], dtype=str))

        conflog_raw_data.append(np.genfromtxt(path + '/' + conflog_files[i], delimiter = ',', usecols = [0, 3, 4, 5, 6, 7]))
        conflog_callsigns.append(np.genfromtxt(path + '/' + conflog_files[i], delimiter = ',', usecols = [1, 2], dtype=str))

    # Determine each single conflict for each vehicle
    t_conflict = [] # Duration of the conflict
    min_dist_conflict = [] # Minimal dist encountered during conflict
    min_t_los_conflict = [] # Minimal time to los during conflict
    for i in range(len(conflog_raw_data)):
        # Determine rows for which UAv0 is the ownship
        conflict_data = filter_raw_conflict_data(conflog_raw_data[i], conflog_callsigns[i], 'UAV0')

        # Analyse each single conflict occured in scenario
        t_conflict.append([])
        min_dist_conflict.append([])
        min_t_los_conflict.append([])
        for j in range(len(conflict_data)):
            t_conflict[i].append(len(conflict_data[j]) * log_interval)

            dist = []
            t_los = []
            for k in range(len(conflict_data[j])):
                dist.append(conflict_data[j][k][2])
                t_los.append(conflict_data[j][k][5])
            min_dist_conflict[i].append(min(dist))
            min_t_los_conflict[i].append(min(t_los))
    return t_conflict, min_dist_conflict, min_t_los_conflict


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

t1, dist1, t_los1 = analyse_ruleset_in_testseries(1, 1, "output")
t3, dist3, t_los3 = analyse_ruleset_in_testseries(1, 3, "output")

flat_dist1 = []
flat_indices1 = []
for i in range(len(dist1)):
    for j in range(len(dist1[i])):
        flat_dist1.append(dist1[i][j])
        flat_indices1.append(i)

flat_dist3 = []
flat_indices3 = []
for i in range(len(dist3)):
    for j in range(len(dist3[i])):
        flat_dist3.append(dist3[i][j])
        flat_indices3.append(i)

print(flat_indices1)
plt.scatter(flat_indices3, flat_dist3)
plt.scatter(flat_indices1, flat_dist1)
plt.show()