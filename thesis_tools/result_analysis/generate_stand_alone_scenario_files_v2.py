import numpy as np 
import os
import json
import parameters
import sys
sys.path.insert(1, 'bluesky/tools')
import geo

kts = 0.514444444 # 1 kts = 0.514444 m/s
nm = 1852. # 1 nm = 1852 m

# Check if there is a scenario.json, otherwise stop
if (not os.path.isfile("thesis_tools/data/scenario.json")):
    print("first generate scenarios using scenario_generator.py")
    exit()

# Check if there is a route.json, otherwise stop
if (not os.path.isfile("thesis_tools/data/route.json")):
    print("first generate routes using route_generator.py")
    exit()

# Check if there is a wind.json, otherwise stop
if (not os.path.isfile("thesis_tools/data/wind.json")):
    print("first generate wind using wind_generator.py")
    exit()

# Check if there is a geofence.json, otherwise stop
if (not os.path.isfile("thesis_tools/data/geofence.json")):
    print("first generate geofences using geofence_generator.py")
    exit()

# Read all relevant json files in data folder
# Load scenario.json
scenario_json = open("thesis_tools/data/scenario.json", "r")
scenario_data = json.load(scenario_json)
scenario_json.close()

# Load route.json
route_json = open("thesis_tools/data/route.json", "r")
route_data = json.load(route_json)
route_json.close()

# Load wind.json
wind_json = open("thesis_tools/data/wind.json", "r")
wind_data = json.load(wind_json)
wind_json.close()

# Load geofence.json
geofence_json = open("thesis_tools/data/geofence.json", "r")
geofence_data = json.load(geofence_json)
geofence_json.close()

# Load coordinates of reference ICAO location
icao_data = np.genfromtxt("data/navdata/airports.dat", delimiter=',', dtype=(str))
airport_found = False
lat_ref = 0.
lon_ref = 0.
for i in range(len(icao_data)):
    if (parameters.ref_position == icao_data[i][0]):
        airport_found = True
        lat_ref = float(icao_data[i][2])
        lon_ref = float(icao_data[i][3])
        break
if not airport_found:
    print("reference position not found")
    exit()

# Scenario file locations and create
loc_thesis_scns = "scenario/Thesis/analysis"
loc_TS0 = loc_thesis_scns + "/TS0"
if not os.path.exists(loc_TS0):
    os.makedirs(loc_TS0)
loc_TS1 = loc_thesis_scns + "/TS1"
if not os.path.exists(loc_TS1):
    os.makedirs(loc_TS1)
loc_TS2 = loc_thesis_scns + "/TS2"
if not os.path.exists(loc_TS2):
    os.makedirs(loc_TS2)
loc_TS3 = loc_thesis_scns + "/TS3"
if not os.path.exists(loc_TS3):
    os.makedirs(loc_TS3)
loc_TS4 = loc_thesis_scns + "/TS4"
if not os.path.exists(loc_TS4):
    os.makedirs(loc_TS4)
loc_TS5 = loc_thesis_scns + "/TS5"
if not os.path.exists(loc_TS5):
    os.makedirs(loc_TS5)

# Empty file locations
filelist = [ f for f in os.listdir(loc_TS0) if f.endswith(".scn") ]
for f in filelist:
    os.remove(os.path.join(loc_TS0, f))

filelist = [ f for f in os.listdir(loc_TS1) if f.endswith(".scn") ]
for f in filelist:
    os.remove(os.path.join(loc_TS1, f))

filelist = [ f for f in os.listdir(loc_TS2) if f.endswith(".scn") ]
for f in filelist:
    os.remove(os.path.join(loc_TS2, f))

filelist = [ f for f in os.listdir(loc_TS3) if f.endswith(".scn") ]
for f in filelist:
    os.remove(os.path.join(loc_TS3, f))

filelist = [ f for f in os.listdir(loc_TS4) if f.endswith(".scn") ]
for f in filelist:
    os.remove(os.path.join(loc_TS4, f))
filelist = [ f for f in os.listdir(loc_TS5) if f.endswith(".scn") ]
for f in filelist:
    os.remove(os.path.join(loc_TS5, f))

# Loop over mission
for i in range(parameters.N_missions):
    # Write scenario files
    # Create first UAV in scenario for wind calm test series
    start_lat0_nowind, start_lon0_nowind = geo.qdrpos(lat_ref, lon_ref, scenario_data[i]['start_qdr0_nowind'], scenario_data[i]['start_dist0_nowind'] / nm)
    start_lat1_nowind, start_lon1_nowind = geo.qdrpos(lat_ref, lon_ref, scenario_data[i]['start_qdr1_nowind'], scenario_data[i]['start_dist1_nowind'] / nm)
    start_lat0_wind, start_lon0_wind = geo.qdrpos(lat_ref, lon_ref, scenario_data[i]['start_qdr0_wind'], scenario_data[i]['start_dist0_wind'] / nm)
    start_lat1_wind, start_lon1_wind = geo.qdrpos(lat_ref, lon_ref, scenario_data[i]['start_qdr1_wind'], scenario_data[i]['start_dist1_wind'] / nm)

    cre_line_wind_calm = "00:00:00.00>CRE UAV0 UAV_" + str(i) + "_0 " + str(start_lat0_nowind) + " " + str(start_lon0_nowind) + " " + str(scenario_data[i]['trk0']) + " " +\
            str(parameters.ref_alt) + " " + str(scenario_data[i]['spd0'] / kts) + "\n"

    # Create first UAV in scenario for windy test series
    cre_line_wind = "00:00:00.00>CRE UAV0 UAV_" + str(i) + "_0 " + str(start_lat0_wind) + " " + str(start_lon0_wind) + " " + str(np.rad2deg(scenario_data[i]['hdg0_wind'])) + " " +\
            str(parameters.ref_alt) + " " + str(scenario_data[i]['spd0'] / kts) + "\n"

    # Create conflicitng UAV in scenario for wind calm scenarios
    conf_line_wind_calm = "00:00:00.00>CRE UAV1 UAV_" + str(i) + "_1 " + str(start_lat1_nowind) + " " + str(start_lon1_nowind) + " " + str(scenario_data[i]['trk1']) + " " +\
            str(parameters.ref_alt) + " " + str(scenario_data[i]['spd1'] / kts) + "\n"

    # Create conflicitng UAV in scenario for windy scenarios
    conf_line_wind = "00:00:00.00>CRE UAV1 UAV_" + str(i) + "_1 " + str(start_lat1_wind) + " " + str(start_lon1_wind) + " " + str(np.rad2deg(scenario_data[i]['hdg1_wind'])) + " " +\
            str(parameters.ref_alt) + " " + str(scenario_data[i]['spd1'] / kts) + "\n"

    # Set bank limits for UAVs
    bank_limit_lines =  "00:00:00.00>BANK UAV0 " + str(parameters.max_bank_angle) + "\n" +\
                        "00:00:00.00>BANK UAV1 " + str(parameters.max_bank_angle) + "\n"

    flyturn_lines = "00:00:00.00>ADDWPT UAV0 FLYOVER\n" +\
                    "00:00:00.00>ADDWPT UAV1 FLYOVER\n"

    # Create routes
    wpt_lines = ""
    for j in range(0, len(route_data[i][0]['points'])):
        wp_lat, wp_lon = geo.qdrpos(lat_ref, lon_ref, route_data[i][0]['qdr'][j], route_data[i][0]['dist'][j] / nm)
        wpt_lines += "00:00:00.00>ADDWPT UAV0 " + str(wp_lat) + "," + str(wp_lon) + "\n"

    for j in range(0, len(route_data[i][1]['points'])):
        wp_lat, wp_lon = geo.qdrpos(lat_ref, lon_ref, route_data[i][1]['qdr'][j], route_data[i][1]['dist'][j] / nm)
        wpt_lines += "00:00:00.00>ADDWPT UAV1 " + str(wp_lat) + "," + str(wp_lon) + "\n"

    # Create geofence
    gf_lines = "00:00:00.00>POLY GF "
    for j in range(len(geofence_data[i]['points'])):
        gf_lat, gf_lon = geo.qdrpos(lat_ref, lon_ref, geofence_data[i]['qdr'][j], geofence_data[i]['dist'][j] / nm)
        gf_lines += str(gf_lat) + " " + str(gf_lon) + " "
    
    gf_lines += "\n"

    # Create wind
    wind_line = "00:00:00.00>WIND 0 0 1000 " + str(wind_data[i]['direction']) + " " + str(wind_data[i]['speed'] / kts) + "\n"

    for j in range(parameters.N_RS):
        # Create scenario files
        scn_TS0 = open(loc_TS0 + "/test" + "{:06d}".format(i) + "_RS" + str(j+1) + ".scn", "w")
        scn_TS1 = open(loc_TS1 + "/test" + "{:06d}".format(i) + "_RS" + str(j+1) + ".scn", "w")
        scn_TS2 = open(loc_TS2 + "/test" + "{:06d}".format(i) + "_RS" + str(j+1) + ".scn", "w")
        scn_TS3 = open(loc_TS3 + "/test" + "{:06d}".format(i) + "_RS" + str(j+1) + ".scn", "w")
        scn_TS4 = open(loc_TS4 + "/test" + "{:06d}".format(i) + "_RS" + str(j+1) + ".scn", "w")
        scn_TS5 = open(loc_TS5 + "/test" + "{:06d}".format(i) + "_RS" + str(j+1) + ".scn", "w")

        common_lines = ""
        common_lines += "00:00:00.00>ASAS ON\n"
        common_lines += "00:00:00.00>RESO SSDUAV\n"
        common_lines += "00:00:00.00>PRIORULES ON RS" + str(j + 1) + "\n"

        # Write scenario files
        scn_TS0.write(cre_line_wind_calm + conf_line_wind_calm + bank_limit_lines + flyturn_lines + wpt_lines + gf_lines)
        scn_TS1.write(common_lines + cre_line_wind_calm + conf_line_wind_calm + bank_limit_lines + flyturn_lines + wpt_lines)
        scn_TS2.write(common_lines + cre_line_wind_calm + conf_line_wind_calm + bank_limit_lines + flyturn_lines + wpt_lines + gf_lines)
        scn_TS3.write(common_lines + wind_line + cre_line_wind + conf_line_wind + bank_limit_lines + flyturn_lines + wpt_lines)
        scn_TS4.write(common_lines + wind_line + cre_line_wind + conf_line_wind + bank_limit_lines + flyturn_lines + wpt_lines + gf_lines)
        scn_TS5.write(wind_line + cre_line_wind + conf_line_wind + bank_limit_lines + flyturn_lines + wpt_lines + gf_lines)

        # Close scenario files
        scn_TS0.close()
        scn_TS1.close()
        scn_TS2.close()
        scn_TS3.close()
        scn_TS4.close()
        scn_TS5.close()