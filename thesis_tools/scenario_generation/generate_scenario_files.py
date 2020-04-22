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

# Check if there is a route_wind.json, otherwise stop
if (not os.path.isfile("thesis_tools/data/route_wind.json")):
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

# Check if there is a geofence_wind.json, otherwise stop
if (not os.path.isfile("thesis_tools/data/geofence_wind.json")):
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

# Load route.json
route_wind_json = open("thesis_tools/data/route_wind.json", "r")
route_wind_data = json.load(route_wind_json)
route_wind_json.close()

# Load wind.json
wind_json = open("thesis_tools/data/wind.json", "r")
wind_data = json.load(wind_json)
wind_json.close()

# Load geofence.json
geofence_json = open("thesis_tools/data/geofence.json", "r")
geofence_data = json.load(geofence_json)
geofence_json.close()

# Load geofence.json
geofence_wind_json = open("thesis_tools/data/geofence_wind.json", "r")
geofence_wind_data = json.load(geofence_wind_json)
geofence_wind_json.close()

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
loc_thesis_scns = "scenario/Thesis"
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

# Empty file locations
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

# Create batch files
batch_TS1 = open(loc_TS1 + "/batch.scn", "w")
batch_TS2 = open(loc_TS2 + "/batch.scn", "w")
batch_TS3 = open(loc_TS3 + "/batch.scn", "w")
batch_TS4 = open(loc_TS4 + "/batch.scn", "w")

# Loop over mission
for i in range(parameters.N_missions):
    # Create scenario files
    scn_TS1 = open(loc_TS1 + "/test" + str(i) + ".scn", "w")
    scn_TS2 = open(loc_TS2 + "/test" + str(i) + ".scn", "w")
    scn_TS3 = open(loc_TS3 + "/test" + str(i) + ".scn", "w")
    scn_TS4 = open(loc_TS4 + "/test" + str(i) + ".scn", "w")

    # Write scenario files
    # Create first UAV in scenario for wind calm test series
    cre_line_wind_calm = "00:00:00.00>CRE UAV0 UAV_" + str(i) + "_0 " + parameters.ref_position + " " + str(scenario_data[i]['trk0']) + " " +\
            str(parameters.ref_alt) + " " + str(scenario_data[i]['spd0'] / kts) + "\n"

    # Create first UAV in scenario for windy test series
    cre_line_wind = "00:00:00.00>CRE UAV0 UAV_" + str(i) + "_0 " + parameters.ref_position + " " + str(scenario_data[i]['hdg0_wind']) + " " +\
            str(parameters.ref_alt) + " " + str(scenario_data[i]['spd0'] / kts) + "\n"
    
    # Create conflicitng UAV in scenario for wind calm scenarios
    conf_line_wind_calm = "00:00:00.00>CRECONFS UAV1 UAV_" + str(i) + "_1 UAV0 " + str(scenario_data[i]['d_psi']) + " " + str(scenario_data[i]['dist_cpa'] / nm) + " " +\
            str(parameters.t_la + parameters.t_extra) + " 0 0 " + str(scenario_data[i]['spd1'] / kts) + "\n"

    # Create conflicitng UAV in scenario for windy scenarios
    conf_line_wind = "00:00:00.00>CRECONFS UAV1 UAV_" + str(i) + "_1 UAV0 " + str(scenario_data[i]['d_psi']) + " " + str(scenario_data[i]['dist_cpa'] / nm) + " " +\
            str(parameters.t_la + parameters.t_extra) + " 0 0 " + str(route_wind_data[i][1]['gs'][0] / kts) + "\n"

    flyturn_line = "00:00:00.00>ADDWPT UAV0 FLYOVER\n"

    # Create routes for wind calm scenarios
    wpt_lines_wind_calm = ""
    for j in range(1, len(route_data[i][0]['points'])):
        wp_lat, wp_lon = geo.qdrpos(lat_ref, lon_ref, route_data[i][0]['qdr'][j], route_data[i][0]['dist'][j] / nm)
        wpt_lines_wind_calm += "00:00:00.00>ADDWPT UAV0 " + str(wp_lat) + "," + str(wp_lon) + "\n"

    for j in range(1, len(route_data[i][1]['points'])):
        wp_lat, wp_lon = geo.qdrpos(lat_ref, lon_ref, route_data[i][1]['qdr'][j], route_data[i][1]['dist'][j] / nm)
        wpt_lines_wind_calm += "00:00:00.00>ADDWPT UAV1 " + str(wp_lat) + "," + str(wp_lon) + "\n"

    # Create routes for windy scenarios
    wpt_lines_wind = ""
    for j in range(1, len(route_wind_data[i][0]['points'])):
        wp_lat, wp_lon = geo.qdrpos(lat_ref, lon_ref, route_wind_data[i][0]['qdr'][j], route_wind_data[i][0]['dist'][j] / nm)
        wpt_lines_wind += "00:00:00.00>ADDWPT UAV0 " + str(wp_lat) + "," + str(wp_lon) + "\n"

    for j in range(1, len(route_data[i][1]['points'])):
        wp_lat, wp_lon = geo.qdrpos(lat_ref, lon_ref, route_wind_data[i][1]['qdr'][j], route_wind_data[i][1]['dist'][j] / nm)
        wpt_lines_wind += "00:00:00.00>ADDWPT UAV1 " + str(wp_lat) + "," + str(wp_lon) + "\n"

    # Create geofences for wind calm scenarios
    gf_lines_wind_calm = "00:00:00.00>POLY GF "
    for j in range(len(geofence_data[i]['points'])):
        gf_lat, gf_lon = geo.qdrpos(lat_ref, lon_ref, geofence_data[i]['qdr'][j], geofence_data[i]['dist'][j] / nm)
        gf_lines_wind_calm += str(gf_lat) + " " + str(gf_lon) + " "
    
    gf_lines_wind_calm += "\n"

    # Create geofences for windy scenarios
    gf_lines_wind = "00:00:00.00>POLY GF "
    for j in range(len(geofence_wind_data[i]['points'])):
        gf_lat, gf_lon = geo.qdrpos(lat_ref, lon_ref, geofence_wind_data[i]['qdr'][j], geofence_wind_data[i]['dist'][j] / nm)
        gf_lines_wind += str(gf_lat) + " " + str(gf_lon) + " "
    
    gf_lines_wind += "\n"

    # Create wind
    wind_line = "00:00:00.00>WIND 0 0 1000 " + str(wind_data[i]['direction']) + " " + str(wind_data[i]['speed'] / kts) + "\n"

    # Write scenario files
    scn_TS1.write(cre_line_wind_calm + conf_line_wind_calm + flyturn_line + wpt_lines_wind_calm)
    scn_TS2.write(cre_line_wind_calm + conf_line_wind_calm + flyturn_line + wpt_lines_wind_calm + gf_lines_wind_calm)
    scn_TS3.write(wind_line + cre_line_wind + conf_line_wind + flyturn_line + wpt_lines_wind)
    scn_TS4.write(wind_line + cre_line_wind + conf_line_wind + flyturn_line + wpt_lines_wind + gf_lines_wind)

    for j in range(parameters.N_RS):
        batch_TS1.write("00:00:00.00>SCEN test_" + str(i) + "_TS1_RS" + str(j + 1) + "\n")
        batch_TS2.write("00:00:00.00>SCEN test_" + str(i) + "_TS2_RS" + str(j + 1) + "\n")
        batch_TS3.write("00:00:00.00>SCEN test_" + str(i) + "_TS3_RS" + str(j + 1) + "\n")
        batch_TS4.write("00:00:00.00>SCEN test_" + str(i) + "_TS4_RS" + str(j + 1) + "\n")

        common_lines = ""
        common_lines += "00:00:00.00>ASAS ON\n"
        common_lines += "00:00:00.00>RESO SSDUAV\n"
        common_lines += "00:00:00.00>PRIORULES ON RS" + str(j + 1) + "\n"

        batch_TS1.write(common_lines)
        batch_TS2.write(common_lines)
        batch_TS3.write(common_lines)
        batch_TS4.write(common_lines)

        batch_TS1.write("00:00:00.00>PCALL Thesis/TS1/test" + str(i) + ".scn\n")
        batch_TS2.write("00:00:00.00>PCALL Thesis/TS2/test" + str(i) + ".scn\n")
        batch_TS3.write("00:00:00.00>PCALL Thesis/TS3/test" + str(i) + ".scn\n")
        batch_TS4.write("00:00:00.00>PCALL Thesis/TS4/test" + str(i) + ".scn\n")

        # Use logging plugin to log
        batch_TS1.write("00:00:00.00>INIT_LOGGERS 1 " + str(j + 1) + "\n")
        batch_TS2.write("00:00:00.00>INIT_LOGGERS 2 " + str(j + 1) + "\n")
        batch_TS3.write("00:00:00.00>INIT_LOGGERS 3 " + str(j + 1) + "\n")
        batch_TS4.write("00:00:00.00>INIT_LOGGERS 4 " + str(j + 1) + "\n")

        # Fast forward
        batch_TS1.write("00:00:00.00>FF\n")
        batch_TS2.write("00:00:00.00>FF\n")
        batch_TS3.write("00:00:00.00>FF\n")
        batch_TS4.write("00:00:00.00>FF\n")

        # end the loggers and save file
        batch_TS1.write("00:00:00.00>SCHEDULE 00:03:00.00 STOP_LOGGERS\n")
        batch_TS2.write("00:00:00.00>SCHEDULE 00:03:00.00 STOP_LOGGERS\n")
        batch_TS3.write("00:00:00.00>SCHEDULE 00:03:00.00 STOP_LOGGERS\n")
        batch_TS4.write("00:00:00.00>SCHEDULE 00:03:00.00 STOP_LOGGERS\n")

        batch_TS1.write("00:00:00.00>SCHEDULE 00:03:00.00 HOLD\n")
        batch_TS2.write("00:00:00.00>SCHEDULE 00:03:00.00 HOLD\n")
        batch_TS3.write("00:00:00.00>SCHEDULE 00:03:00.00 HOLD\n")
        batch_TS4.write("00:00:00.00>SCHEDULE 00:03:00.00 HOLD\n")

    # Close scenario files
    scn_TS1.close()
    scn_TS2.close()
    scn_TS3.close()
    scn_TS4.close()

# Close batch files
batch_TS1.close()
batch_TS2.close()
batch_TS3.close()
batch_TS4.close()
